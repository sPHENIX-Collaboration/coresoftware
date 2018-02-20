```g4cemc``` is deprecated. Components of this lib has been distributed to various new locations following the following map. Further discussion see https://github.com/sPHENIX-Collaboration/coresoftware/pull/403

![relocation map](https://user-images.githubusercontent.com/7947083/34732682-f3cf99c2-f533-11e7-883c-f945d9871ddf.png)

<details> 
<summary></summary>
relocation_map
digraph G {
    ordering="in";
    rankdir="LR";
    splines="line";
    dpi="300";
    edge [comment="Wildcard node added automatic in EG."];
    node [comment="Wildcard node added automatic in EG."];
    subgraph cluster_g4cemc {
        bgcolor="coral";
        label="simulation/g4simulation/g4cemc";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        g4cemc_RawTowerDefs [label="RawTowerDefs.h"];
        g4cemc_RawTower [label="RawTower*"];
        g4cemc_RawTowerGeom [label="RawTowerGeom*"];
        g4cemc_RawCluster [label="RawCluster*"];
        g4cemc_RawTowerBuilder [label="RawTowerBuilder*"];
        g4cemc_RawTowerDigitizer [label="RawTowerDigitizer*"];
        g4cemc_RawClusterBuilder [color="blue4", 
                                  fontcolor="blue4", 
                                  shape="box", 
                                  fillcolor="white", 
                                  label="RawClusterBuilder"];
        g4cemc_RawClusterBuilderv1 [color="blue4", 
                                    fontcolor="blue4", 
                                    shape="box", 
                                    fillcolor="white", 
                                    label="RawClusterBuilderv1"];
        g4cemc_RawClusterPositionCorrection [label="RawClusterPositionCorrection*"];
        g4cemc_RawTowerCalibration [label="RawTowerCalibration*"];
        g4cemc_RawTowerCombiner [label="RawTowerCombiner*"];
    }

    subgraph cluster_g4calo {
        bgcolor="AliceBlue";
        label="simulation/g4simulation/g4calo";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        g4calo_RawTowerBuilder [label="RawTowerBuilder*"];
        g4calo_RawTowerDigitizer [label="RawTowerDigitizer*"];
        libg4calo [bgcolor="coral", 
                   shape="box", 
                   fillcolor="white", 
                   color="aquamarine4", 
                   fontcolor="aquamarine4", 
                   label="libg4calo.so"];
        g4calo_RawTowerBuilder -> libg4calo  [color="aquamarine4"];
        g4calo_RawTowerDigitizer -> libg4calo  [color="aquamarine4"];
    }

    subgraph cluster_calobase {
        bgcolor="AliceBlue";
        label="offline/package/CaloBase";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        calobase_RawTowerDefs [label="RawTowerDefs.h"];
        calobase_RawTower [label="RawTower*"];
        calobase_RawTowerGeom [label="RawTowerGeom*"];
        calobase_RawCluster [label="RawCluster*"];
        calobase_RawClusterUtility [color="aquamarine4", 
                                    fontcolor="aquamarine4", 
                                    label="RawClusterUtility"];
        libcalo_io [bgcolor="coral", 
                    shape="box", 
                    fillcolor="white", 
                    color="aquamarine4", 
                    fontcolor="aquamarine4", 
                    label="libcalo_io.so"];
        libcalo_util [bgcolor="coral", 
                      shape="box", 
                      fillcolor="white", 
                      color="aquamarine4", 
                      fontcolor="aquamarine4", 
                      label="libcalo_util.so"];
        calobase_RawTowerDefs -> libcalo_io  [color="aquamarine4"];
        calobase_RawTower -> libcalo_io  [color="aquamarine4"];
        calobase_RawTowerGeom -> libcalo_io  [color="aquamarine4"];
        calobase_RawCluster -> libcalo_io  [color="aquamarine4"];
        calobase_RawClusterUtility -> libcalo_util  [color="aquamarine4"];
    }

    subgraph cluster_caloreco {
        bgcolor="AliceBlue";
        label="offline/package/CaloReco";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        caloreco_RawClusterBuilderGraph [color="blue4", 
                                         fontcolor="blue4", 
                                         shape="box", 
                                         fillcolor="white", 
                                         label="RawClusterBuilderGraph"];
        caloreco_RawClusterBuilderTemplate [color="blue4", 
                                            fontcolor="blue4", 
                                            shape="box", 
                                            fillcolor="white", 
                                            label="RawClusterBuilderTemplate"];
        caloreco_RawClusterPositionCorrection [label="RawClusterPositionCorrection*"];
        caloreco_RawTowerCalibration [label="RawTowerCalibration*"];
        caloreco_RawTowerCombiner [label="RawTowerCombiner*"];
        libcalo_reco [bgcolor="coral", 
                      shape="box", 
                      fillcolor="white", 
                      color="aquamarine4", 
                      fontcolor="aquamarine4", 
                      label="libcalo_reco.so"];
        caloreco_RawClusterBuilderGraph -> libcalo_reco  [color="aquamarine4"];
        caloreco_RawClusterBuilderTemplate -> libcalo_reco  [color="aquamarine4"];
        caloreco_RawClusterPositionCorrection -> libcalo_reco  [color="aquamarine4"];
        caloreco_RawTowerCalibration -> libcalo_reco  [color="aquamarine4"];
        caloreco_RawTowerCombiner -> libcalo_reco  [color="aquamarine4"];
    }

    subgraph cluster_legend {
        color="white";
        bgcolor="gray94";
        label="Legend";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        Moved [label="Moved"];
        Renamed [style="filled", 
                 fillcolor="white", 
                 color="blue4", 
                 fontcolor="blue4", 
                 shape="box", 
                 label="Renamed"];
        Removed [bgcolor="coral", 
                 shape="box", 
                 label="Removed ", 
                 style="filled", 
                 fillcolor="coral"];
        New [bgcolor="coral", 
             shape="box", 
             fillcolor="white", 
             color="aquamarine4", 
             fontcolor="aquamarine4", 
             label="New"];
    }

    g4cemc_RawTowerDefs -> calobase_RawTowerDefs;
    g4cemc_RawTower -> calobase_RawTower;
    g4cemc_RawTowerGeom -> calobase_RawTowerGeom;
    g4cemc_RawCluster -> calobase_RawCluster;
    g4cemc_RawTowerBuilder -> g4calo_RawTowerBuilder;
    g4cemc_RawTowerDigitizer -> g4calo_RawTowerDigitizer;
    g4cemc_RawClusterBuilder -> caloreco_RawClusterBuilderGraph  [color="blue4"];
    g4cemc_RawClusterBuilderv1 -> caloreco_RawClusterBuilderTemplate  [color="blue4"];
    g4cemc_RawClusterPositionCorrection -> caloreco_RawClusterPositionCorrection;
    g4cemc_RawTowerCalibration -> caloreco_RawTowerCalibration;
    g4cemc_RawTowerCombiner -> caloreco_RawTowerCombiner;
}
relocation_map
</details>

![Alt text](https://g.gravizo.com/source/relocation_map?https%3A%2F%2Fraw.githubusercontent.com%2Fblackcathj%2Fcoresoftware%2FClusterizerRelocated%2Fsimulation%2Fg4simulation%2Fg4cemc%2FREADME.md)
