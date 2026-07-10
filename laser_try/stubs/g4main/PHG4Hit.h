#pragma once
// host stub: minimal PHG4Hit base for running collaboration generators verbatim
#ifndef STUB_NO_UNITS
static constexpr double mm = 1.0;
static constexpr double cm = 10.0;
#endif
class PHG4Hit
{
 public:
  virtual ~PHG4Hit() = default;
  double x[2] = {0, 0}, y[2] = {0, 0}, z[2] = {0, 0}, t[2] = {0, 0};
  double px[2] = {0, 0}, py[2] = {0, 0}, pz[2] = {0, 0};
  double edep = 0, eion = 0;
  int layer = 0, trkid = 0;
  void set_x(int i, double v) { x[i] = v; }
  void set_y(int i, double v) { y[i] = v; }
  void set_z(int i, double v) { z[i] = v; }
  void set_t(int i, double v) { t[i] = v; }
  void set_px(int i, double v) { px[i] = v; }
  void set_py(int i, double v) { py[i] = v; }
  void set_pz(int i, double v) { pz[i] = v; }
  void set_layer(int v) { layer = v; }
  void set_trkid(int v) { trkid = v; }
  void set_edep(double v) { edep = v; }
  void set_eion(double v) { eion = v; }
  double get_x(int i) const { return x[i]; }
  double get_y(int i) const { return y[i]; }
  double get_z(int i) const { return z[i]; }
  double get_t(int i) const { return t[i]; }
  double get_edep() const { return edep; }
  double get_eion() const { return eion; }
  int get_hit_id() const { return 0; }
};
