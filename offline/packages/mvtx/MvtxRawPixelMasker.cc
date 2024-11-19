#include "MvtxRawPixelMasker.h"
#include "MvtxPixelMask.h"
#include "MvtxPixelMaskv1.h"
#include "MvtxRawPixelDefs.h"

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHTypedNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <ffarawobjects/MvtxRawEvtHeader.h>
#include <ffarawobjects/MvtxRawEvtHeaderv2.h>
#include <ffarawobjects/MvtxRawHit.h>
#include <ffarawobjects/MvtxRawHitv1.h>
#include <ffarawobjects/MvtxRawHitContainer.h>
#include <ffarawobjects/MvtxRawHitContainerv1.h>

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdlib>

int MvtxRawPixelMasker::InitRun(PHCompositeNode* topNode)
{
    std::cout << "MvtxRawPixelMasker::InitRun - Initializing" << std::endl;
    m_current_pixel_hits.clear();
    return CreateNode(topNode);
}


int MvtxRawPixelMasker::process_event(PHCompositeNode *topNode)
{
    if( !m_dynamic_mode ) { return Fun4AllReturnCodes::EVENT_OK; }
    
    // get dst nodes
    auto mvtx_raw_event_header = findNode::getClass<MvtxRawEvtHeaderv2>(topNode, "MVTXRAWEVTHEADER");
    if(!mvtx_raw_event_header){
        std::cout << PHWHERE << "::" << __func__ << ": Could not get MVTXRAWEVTHEADER from Node Tree" << std::endl;
        exit(1); 
    }
    if (Verbosity() > 2) {  mvtx_raw_event_header->identify();  }

    auto mvtx_raw_hit_container = findNode::getClass<MvtxRawHitContainerv1>(topNode, "MVTXRAWHIT");
    if (!mvtx_raw_hit_container){ 
        std::cout << PHWHERE << "::" << __func__ << ": Could not get MVTXRAWHIT from Node Tree" << std::endl; 
        exit(1); 
    }
    if(Verbosity() > 2) {  mvtx_raw_hit_container->identify();  }

    auto mvt = findNode::getClass<MvtxPixelMaskv1>(topNode, m_hot_pixel_node_name);
    if (!mvt){
        std::cout << PHWHERE << "::" << __func__ << ": Could not get MvtxPixelMaskv1 from Node Tree" << std::endl;
        exit(1);
    }
    if (Verbosity() > 2) { mvt->identify(); }

    // get LL1s from the MVTX raw event header
    for (unsigned int ihit = 0; ihit < mvtx_raw_hit_container->get_nhits(); ihit++)
    {
        // get this hit
        auto mvtx_hit = mvtx_raw_hit_container->get_hit(ihit);
        if (!mvtx_hit){ 
            std::cout << PHWHERE << "::" << __func__ << ": Could not get MVTX hit from container. Hit index: " << ihit << std::endl;
            continue; 
        }
        if ( mvt->is_masked(mvtx_hit) ){ // skip previously masked pixels
            continue;
        }

        // get the hit info
        uint64_t strobe = mvtx_hit->get_bco();

        if ( strobe > m_last_strobe ){ 
            m_last_strobe = strobe; 
            m_current_strobe_counter++;
        }

        addHit(mvtx_hit);        
    }

    if (m_current_strobe_counter > m_strobe_update)
    {
        if(Verbosity() > 2) { 
            std::cout << "MvtxRawPixelMasker::process_event - Updating mask" << std::endl; 
        }
        // Sort the pixel hit vector
        std::sort(m_current_pixel_hits.begin(), m_current_pixel_hits.end(), [](const MvtxRawPixelDefs::pixelhitpair_t& a, const MvtxRawPixelDefs::pixelhitpair_t& b) { return a.second > b.second; });

        // get initial values
        unsigned int nhits_all = 0;
        for(auto it = m_current_pixel_hits.begin(); it != m_current_pixel_hits.end(); ++it) {
            nhits_all += it->second;
        }
        unsigned int nHits_threshold = static_cast<unsigned int>(1.0*MvtxRawPixelDefs::kNpixelsTotal*m_current_strobe_counter*m_threshold);
        
        int ipixel = 0;
        while(nhits_all > nHits_threshold)
        {
            unsigned int nhits = m_current_pixel_hits[ipixel].second;
            MvtxRawPixelDefs::mvtx_pixelkey key = m_current_pixel_hits[ipixel].first;
            nhits_all -= nhits;
            mvt->add_pixel(key);
            m_npixels_masked++;
            ipixel++;
        }

        m_current_strobe_counter = 0;
        m_current_pixel_hits.clear();
    }

    return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxRawPixelMasker::End(PHCompositeNode * /*topNode*/)
{
    std::cout << "MvtxRawPixelMasker::End - Masked " << m_npixels_masked << " pixels" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;

}

int MvtxRawPixelMasker::CreateNode(PHCompositeNode *topNode)
{
    PHNodeIterator iter(topNode);
    auto runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
    if ( !runNode ) {
        std::cout << PHWHERE << "RUN Node missing, doing nothing." << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
    }

    PHNodeIterator iter_run(runNode);
    auto mvtxNode = dynamic_cast<PHCompositeNode *>(iter_run.findFirst("PHCompositeNode", "MVTX"));
    if ( !mvtxNode )
    {
        mvtxNode = new PHCompositeNode("MVTX");
        runNode->addNode(mvtxNode);
    }

    // Create the Input node if required
    MvtxPixelMaskv1 *mvt = findNode::getClass<MvtxPixelMaskv1>(runNode, m_hot_pixel_node_name);
    if (!mvt)
    {
        mvt = new MvtxPixelMaskv1();
        PHIODataNode<PHObject> *MvtxPixelMaskNode = new PHIODataNode<PHObject>(mvt, m_hot_pixel_node_name, "PHObject");
        mvtxNode->addNode(MvtxPixelMaskNode);
    }
    if (m_load_from_cdb) { mvt->load_from_CDB(m_calibfile); }
    auto currentmask = mvt->get_pixel_mask();
    m_npixels_masked = currentmask.size();
    std::cout << "MvtxRawPixelMasker::CreateNode - Created MvtxPixelMaskv1 node" << std::endl;

    return Fun4AllReturnCodes::EVENT_OK;
}

void MvtxRawPixelMasker::addHit(MvtxRawHit *hit, const unsigned int nhits)
{
    // Add a hit to the hit map
    MvtxRawPixelDefs::mvtx_pixelkey key = MvtxRawPixelDefs::gen_pixelkey(hit);
    // see if this pixel is already in the pixel hit pair vector
    auto it = std::find_if(m_current_pixel_hits.begin(), m_current_pixel_hits.end(), [key](const MvtxRawPixelDefs::pixelhitpair_t& pair) { return pair.first == key; });
    if(it != m_current_pixel_hits.end()) {
        it->second+=nhits;
    } else{
        m_current_pixel_hits.push_back(std::make_pair(key, nhits));
    }

    return;
}