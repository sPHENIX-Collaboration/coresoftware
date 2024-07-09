#include <Rtypes.h>
#include <TAxis.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TPaveStats.h>

#include <iostream>
#include <utility>

namespace InttQa
{
  /*!
    @brief A common header file for INTT QA
   */
  // common variables
  const int kFelix_num = 8;     //! the number of our FELIX server
  const int kFee_num = 14;      //! the number of half-ladders in a single FELIX server
  const int kChip_num = 26;     //! the number of chip in a half-ladder
  const int kChan_num = 128;    //! the number of channel in a single chip
  const int kFirst_pid = 3001;  //! the first pid (packet ID), which means intt0
  const int kBco_max = 128;

  const int kColors[10] = {
      kBlack, kRed, kBlue,
      kGreen + 2, kMagenta + 1, kYellow + 1,
      kCyan + 1, kOrange + 1, kBlue + 9,
      kGray + 2};  //! A list of nice colors

  // functions depending on variables above

  //! It returns int, which can be used as color in ROOT, for num-th graph/histograms. It checks if given parameter is in the range or not.
  template <class aaa>
  int GetColor(aaa num)
  {
    assert(0 <= num && num < 10);
    return InttQa::kColors[num];
  }

  //! Some configuration is done for a better-looking histogram
  template <class TH>
  void HistConfig(TH* hist, int index = 0)
  {
    hist->SetLineColor(InttQa::GetColor(index));
    hist->SetFillColorAlpha(hist->GetLineColor(), 0.1);
  }

  // functions depending on variables/functions above

  /*!
    @brief A wrapper function to draw a palette axis as you like
    @param xmin A relative coordinate of minimum x, which can be from 0 to 1.
    @param xmin A relative coordinate of minimum y, which can be from 0 to 1.
    @param xmax A relative coordinate of maximum x, which can be from 0 to 1.
    @param xmax A relative coordinate of maximum y, which can be from 0 to 1.
   */
  template <typename TH>
  TPaletteAxis* DrawPaletteAxis(TH* hist,
                                double xmin, double ymin,
                                double xmax, double ymax,
                                double label_size = 0.04)

  {
    gPad->Update();
    TPaletteAxis* pal = (TPaletteAxis*) hist->GetListOfFunctions()->FindObject("palette");
    pal->GetAxis()->SetLabelSize(label_size);
    pal->GetAxis()->CenterTitle();

    pal->SetX1NDC(xmin);
    pal->SetX2NDC(xmax);

    pal->SetY1NDC(ymin);
    pal->SetY2NDC(ymax);
    pal->Draw();

    return pal;
  }

  /*!
    @brief A wrapper function to draw a statistical box as you like
    @param hist A pointer of a histogram object, which can be any of TH1, TH2, and TH3.
    @param xmin A relative coordinate of minimum x, which can be from 0 to 1.
    @param xmin A relative coordinate of minimum y, which can be from 0 to 1.
    @param xmax A relative coordinate of maximum x, which can be from 0 to 1.
    @param xmax A relative coordinate of maximum y, which can be from 0 to 1.
   */
  template <typename TH>
  void DrawStats(TH* hist, double xmin, double ymin, double xmax, double ymax, int /*font = 4*/)
  {
    gPad->Update();
    TPaveStats* st = (TPaveStats*) hist->FindObject("stats");
    if (st == nullptr)
      return;

    st->SetTextColorAlpha(hist->GetLineColor(), 1.0);
    st->SetLineColorAlpha(hist->GetLineColor(), 1.0);
    st->SetFillStyle(1001);
    st->SetFillColor(0);

    st->SetX1NDC(xmin);
    st->SetX2NDC(xmax);
    st->SetY1NDC(ymin);
    st->SetY2NDC(ymax);

    st->Draw("same");
  }

  //! i don't remember what's this
  template <class TH>
  void HistsConfig(int hist_num, TH* hists)
  {
    std::vector<int> bin_x_with_entry;
    std::vector<int> bin_y_contents;
    for (int i = 0; i < hist_num; i++)
    {
      for (int j = 1; j < hists[i]->GetNbinsX() + 1; j++)
      {
        if (hists[i]->GetBinContent(j) != 0)
        {
          bin_x_with_entry.push_back(j);
          bin_y_contents.push_back(hists[i]->GetBinContent(j));
        }
      }
    }

    int min_non_zero_x = *std::min_element(bin_x_with_entry.begin(), bin_x_with_entry.end());
    if (min_non_zero_x > 1)
      min_non_zero_x--;

    int max_non_zero_x = *std::max_element(bin_x_with_entry.begin(), bin_x_with_entry.end());
    if (max_non_zero_x < hists[0]->GetNbinsX())
      max_non_zero_x++;

    hists[0]->GetXaxis()->SetRange(min_non_zero_x, max_non_zero_x);

    int min_y = *std::min_element(bin_y_contents.begin(), bin_y_contents.end());
    int max_y = *std::max_element(bin_y_contents.begin(), bin_y_contents.end());
    if (gPad->GetLogy() == 0)  // linear
      max_y *= 1.2;
    else
      max_y *= 2;

    hists[0]->GetYaxis()->SetRangeUser(min_y, max_y);
    for (int i = 0; i < hist_num; i++)
    {
      InttQa::HistConfig(hists[i], i);
    }

    std::cout << "X range: " << min_non_zero_x << "\t" << max_non_zero_x << std::endl;
    std::cout << "Y range: " << min_y << "\t" << max_y << std::endl;
  }

  template <class T>
  std::pair<int, int> OptimizeRange(T* hist, int axis_param = 0)
  {
    /*!
      @brief The range of histogram is optimized to show bins with non-zero entries
      @param hist A histogram to be modified
      @param axis_param Choice of axis to be modified. It can be x(0), y(1), or z(2).
      @retval pair < int, int > The minimum bin (.first) and maximum bin (.second)
     */

    TAxis* axis = hist->GetXaxis();
    if (axis_param == 1)
      axis = hist->GetYaxis();
    else if (axis_param == 2)
      axis = hist->GetZaxis();

    int bin_min = 0;                 // container for the minimum bin ID
    int bin_max = axis->GetNbins();  // container for the maximum bin ID
    int bin_max_sweep = bin_max;     // for for loop

    // sweep from low to up to find the minimum bin with non-zero content
    for (int i = 1; i < bin_max_sweep; i++)
    {
      auto content = hist->GetBinContent(i);
      if (content > 0)
      {
        bin_min = i;
        break;
      }
    }

    // sweep from up to low to find the maximum bin with non-zero content
    for (int i = bin_max_sweep - 1; i > 0; i--)
    {
      auto content = hist->GetBinContent(i);
      if (content > 0)
      {
        bin_max = i;
        break;
      }
    }

    axis->SetRange(bin_min, bin_max);
    std::pair<int, int> rtn(bin_min, bin_max);
    return rtn;
  }
}  // end of namespace InttQa
