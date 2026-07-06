// reencode_field.C
// Re-encode a TPC/solenoid field-map "fieldmap" TNtuple in the NATIVE ROOT format of
// the running interpreter. Run this INSIDE the ana.331 container (ROOT 6.24) to
// downgrade a field map written by a far-newer ROOT (e.g. 6.38) that the 6.24 runtime
// cannot stream safely -- reading such a file directly corrupts ROOT's object cleanup
// (TList "already deleted") and SIGSEGVs in the full simulation. Re-filling the ntuple
// entry-by-entry rewrites native 6.24 baskets (a "fast" CloneTree is NOT enough -- it
// copies the original baskets). Values are preserved exactly; only the I/O format changes.
//
// Usage (in container):  root -b -q 'reencode_field.C("<in.root>","<out.root>")'
void reencode_field(const char* in, const char* out)
{
  TFile* fi = TFile::Open(in);
  if (!fi || fi->IsZombie()) { printf("reencode_field ERROR: cannot open %s\n", in); return; }
  TNtuple* t = (TNtuple*) fi->Get("fieldmap");
  if (!t) { printf("reencode_field ERROR: no 'fieldmap' TNtuple in %s\n", in); return; }

  // reproduce the exact variable list (x:y:z:bx:by:bz:hz ...) from the source branches
  TString varlist;
  TIter it(t->GetListOfBranches());
  TBranch* b;
  while ((b = (TBranch*) it())) { if (varlist.Length()) varlist += ":"; varlist += b->GetName(); }

  TFile* fo = TFile::Open(out, "RECREATE");
  TNtuple* o = new TNtuple("fieldmap", "fieldmap", varlist);
  Long64_t n = t->GetEntries();
  for (Long64_t i = 0; i < n; i++) { t->GetEntry(i); o->Fill(t->GetArgs()); }
  o->Write();
  fo->Close();
  fi->Close();
  printf("reencode_field OK: %lld entries  vars=[%s]  -> %s\n", n, varlist.Data(), out);
}
