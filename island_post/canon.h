// canon.h — CANONICAL conventions for every comparison figure/analysis in island_post.
// Conventions live HERE as code, never re-derived from memory in individual macros.
// History of the errors this file exists to prevent: missing layer cut, missing adc>0,
// missing pedestal bridge, half-integer binning on integer ADC, size-semantics mixups.
#pragma once
#include <TFile.h>
#include <TH1D.h>

namespace CANON
{
// Real-side selection for ANY histogram filled from the real ntuplizer file
// (ntp_hit contains ALL trackers; adc=0 = instrumental-burst population, excluded by anchors).
inline const char *TPC_CUT = "layer>=7&&layer<=54&&adc>0";

// Laser-flash fiducial time cut (SHORT-TERM filter, conventional analysis practice).
// Identified as the sPHENIX DIFFUSE LASER calibration flash: 266 nm light liberates
// photoelectrons from the aluminum stripes on the gold central membrane (z=0), so every
// flash arrives ~one full drift later -> fixed tbin ~329, detector-wide, both sides,
// soft ADC, apparent z ~ -35 cm. Collaboration models it: coresoftware
// PHG4TpcCentralMembrane (stripe generator, electrons_per_stripe) + tpccalib
// TpcCentralMembraneMatching (extraction side). 44/100 frames flagged in laser_frames.txt.
// Excluding the fixed window on BOTH datasets identically is unbiased for steady-state
// observables (drops 1.9% of the time axis). LONG-TERM (queued): inject a transported
// PHG4TpcCentralMembrane stripe flash in frame_composer -> labeled noise class cls=2.
inline const char *TPC_CUT_NOLASER = "layer>=7&&layer<=54&&adc>0&&!(tbin>=322&&tbin<=340)";
// (sim digi zbin==tbin exactly, all rows verified 2026-07-20; hit69 exports carry
//  zbin=NaN per the real TPC convention -> cut on tbin, valid for every sim source)
inline const char *SIM_NOLASER = "!(tbin>=322&&tbin<=340)";  // apply symmetrically to sim

// Stock-ana.331 evaluator/port adc and maxadc are pedestal-UNSUBTRACTED; real production
// (and the detached B3 readout) are subtracted. Bridge stock curves with these expressions.
inline const double PEDESTAL = 74.4;
inline const char *STOCK_ADC = "adc-74.4";
inline const char *STOCK_MAXADC = "maxadc-74.4";

// REAL z CONVENTION QUIRK (measured 2026-07-11): real ntp_hit z for zelem==1 is standard
// apparent-z (105.5 - v*t), but for zelem==0 it is v*t (drift distance, offset by +105.5
// vs the physical convention). Convert on ingestion: z_phys = z - 105.5*(zelem==0).
// This artifact faked an 'asymmetric in-time hump' and polluted real apparent-eta.

// Real production ntp_cluster `size` = phisize*zsize (bounding-box AREA; 99.52% exact —
// TrkrClusterv5 keeps no hit lists). Detached islander / port `size` = true pixel count.
// Compare AREA vs AREA (use phisize*zsize on the sim side) or count vs count, never mixed.

// Integer ADC everywhere: histograms of adc MUST use integer-centred bins,
// e.g. (1101,-0.5,1100.5) or (121,-0.5,120.5). Never fractional-width bins.

// Real pixel-level reference curves: CONSUME the frozen anchors; never refill ad hoc.
inline TH1D *anchor(const char *name, const char *ref = "ref_real.root")
{
  TFile *f = TFile::Open(ref);
  if (!f || f->IsZombie())
  {
    printf("CANON::anchor: cannot open %s\n", ref);
    return nullptr;
  }
  TH1D *h = (TH1D *) f->Get(name);
  if (!h)
  {
    printf("CANON::anchor: no %s in %s\n", name, ref);
    return nullptr;
  }
  h = (TH1D *) h->Clone(Form("anchor_%s", name));
  h->SetDirectory(nullptr);
  f->Close();
  return h;
}
}  // namespace CANON
