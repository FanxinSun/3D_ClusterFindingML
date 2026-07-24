void pp_census(const char *fn, const char *chunk, const char *out)
{
  TFile *f = TFile::Open(fn);
  TTree *t = (TTree *) f->Get("ntp_hit");
  float ev; t->SetBranchStatus("*", 0);
  t->SetBranchStatus("event", 1); t->SetBranchAddress("event", &ev);
  std::map<int, long> n;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); n[(int) ev]++; }
  FILE *fo = fopen(out, "a");
  for (auto &kv : n) fprintf(fo, "%s %d %ld\n", chunk, kv.first, kv.second);
  fclose(fo);
}
