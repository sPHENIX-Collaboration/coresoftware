root -l <<EOF
  static const unsigned int kBitShiftClusId = 32;
  uint32_t hskey = 4;
  uint32_t clusid = 15;

  /* pair< uint64_t, pair<uint32_t, uint32_t> > data; */

  class triple {
    public:
    uint32_t hskey;
    uint32_t clusid;
    uint64_t comb;
    triple(uint32_t _hskey, uint32_t _clusid) :
      hskey{_hskey}, clusid{_clusid} {
        uint64_t tmp = _hskey;
        comb = (tmp << kBitShiftClusId);
        comb |= clusid;
      };
  };
  triple a { 23, 23};
  /* bool comp(triple& lhs, triple& rhs) { return lhs.comb < rhs.comb }; */

  vector<triple> dat;
  dat.push_back({ 1, 2 });
  dat.push_back({ 2, 1 });
  dat.push_back({ 0, 10 });
  dat.push_back({ 35, 9 });

  for (auto& D : dat) { cout << setw(12) << D.comb << " " << D.hskey << " " << D.clusid << endl; }
  cout << " sorted: " << endl;
  sort(dat.begin(), dat.end(), [](const triple&lhs, const triple&rhs){return lhs.comb<rhs.comb;});
  for (auto& D : dat) { cout << setw(12) << D.comb << " " << D.hskey << " " << D.clusid << endl; }


  /* const uint64_t tmp = hskey; */
  /* cout << " tmp: " << tmp << endl; */
  /* uint64_t tmp2 = (tmp << kBitShiftClusId); */
  /* cout << " tmp: " << tmp2 << endl; */
  /* tmp2 |= clusid; */
  /* cout << " tmp: " << tmp2 << endl; */

  /* uint32_t _clusid = tmp2; */
  /* cout << " clusid: " << clusid << " vs " << clusid << endl; */
  /* uint32_t _hskey = (tmp2 >> kBitShiftClusId); */
  /* cout << " hskey:  " << hskey << " vs " << _hskey << endl; */
EOF
