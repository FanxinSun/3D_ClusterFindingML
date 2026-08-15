# p5_display — complete option and command reference

Framing note: `/p5/view` and `/p5/zoom` are expressed as "show this radius in
cm" and derive the Geant4 zoom factor from the live scene extent. A raw
`/vis/viewer/zoomTo <f>` is relative to the scene's bounding sphere, which Geant4
computes from the VISIBLE volumes — so a fixed factor frames the TPC under
preset `tpc` and the whole beam line under `all`.

Generated from the source, every entry exercised once in a batch run
(`--vis TSG_OFFSCREEN`) with zero command errors. `p5_display.sh` passes any
unrecognised flag straight through to the binary, so the two lists below are the
whole surface.

---

## 1. Shell options — `./p5_display.sh <flags>`

Run bare (`./p5_display.sh`) for the interactive menu instead.

### Selection

| flag | default | meaning |
|---|---|---|
| `--ver <tag>` | `v55` | production tag; picks `digi_frames_production_<tag>.root` and `island91_frames_production_<tag>.root` |
| `--frame <n>` | `0` | production frame (0..249 for the 5x50 layout) |
| `--g4event <n>` | `-1` = off | show a **library** G4 collision instead of a frame |
| `--g4chunk <i>` | `0` | which `PP_g4hit_<i>.root` `--g4event` reads |

### Filters at startup

| flag | default | meaning |
|---|---|---|
| `--sector <0..11>` | `-1` = all 12 | one TPC 30-degree wedge (~1/12 the content) |
| `--layers <lo> <hi>` | `7 54` | TPC row range |
| `--adcmin <v>` | `0` | drop pixels below this ADC |
| `--pixsize <px>` | `0` frame / `2.5` library | `0` = 1-px dots, `>0` = circles |
| `--clussize <px>` | `2.2` | base screen size of a cluster mark |
| `--zfull` | off | start with the full +-310 cm apparent-z window instead of the physical TPC |

### Look

| flag | default | meaning |
|---|---|---|
| `--colour <mode>` | `level` | `level` (hit vs cluster vs chain) / `class` (ML truth label) / `track` (spotlight) |
| `--geometry <p>` | `tpc` | `tpc` / `tracker` / `calo` / `all` / `none` |
| `--geomstyle <s> [a]` | `auto`, alpha 1.0 | `wireframe` / `surface [alpha]` (alpha defaults 0.12) / `auto` |
| `--view <preset>` | `3d` | `3d` / `rphi` / `rz` / `iso` |
| `--zoom <radius_cm>` | per preset | frame this radius instead of the preset's |
| `--geombright <f>` | `1.0` | geometry brightness multiplier |
| `--geomtint <bool>` | `true` | colour the geometry per detector system (false = flat grey) |

### Paths (all derived from `--repo` + `--ver` unless overridden)

| flag | default | meaning |
|---|---|---|
| `--repo <dir>` | `/home/rog/sPHENIX/3D_ClusterFindingML` | repository root |
| `--digi <file>` | derived | pixel source (`ntp_hit`) |
| `--island <file>` | derived | cluster source (`ntp_cluster` + `ntp_truth`) |
| `--gdml <file>` | `P5/sphenix_p5.gdml` | geometry |
| `--hepmc <file>` | none (particle gun) | HepMC2 chunk feeding `/run/beamOn` |
| `--field <file\|off>` | off | field map; lazy unless `P5_FIELD_PRELOAD` is set |

### Session / output

| flag | default | meaning |
|---|---|---|
| `--vis <driver>` | `OGL` (= `TSG_QT_GLES` in G4 11.4) | also `OGLSQt`, `OGLSX`, `TSG_OFFSCREEN` |
| `--session <type>` | `qt` | `qt` / `tcsh` / `csh` |
| `--size <WxH>` | `1400x900` | window size hint |
| `--macro <file>` | none | extra macro executed after the vis setup |
| `--batch` | off | run and exit, no interactive session |
| `--cmd "<command>"` | — | **any** `/p5/` or `/vis/` command, repeatable, applied in order once the viewer exists |
| `--help` | — | this flag list |

### Environment

| variable | effect |
|---|---|
| `P5_FIELD_PRELOAD` | load the field map at `Initialize()` instead of lazily |
| `P5_DISP_CACHE` | directory for the per-event range index (default `<datadir>/.p5disp`) |

`--cmd` is the general escape hatch: anything in section 2 that has no dedicated
flag is reachable as `--cmd "/p5/whatever ..."`.

---

## 2. GUI commands — 34 commands in 4 directories

Type `help /p5/` in the session, or `/control/manual /p5/` for the full dump
with parameter types.

### Selection

| command | default | meaning |
|---|---|---|
| `/p5/frame <n>` | `0` | load a production frame (also leaves library mode) |
| `/p5/g4event <n>` | `-1` | load a library G4 event; `-1` returns to frame mode |
| `/p5/g4chunk <i>` | `0` | `PP_g4hit` chunk for `/p5/g4event` |
| `/p5/g4list <n>` | — | the N busiest library events in the current chunk, with the multiple of the median |
| `/p5/reload` | — | re-read the current selection |
| `/p5/print` | — | full selection summary (sources, cuts, counts) |

### Filters — these re-read the files

| command | default | meaning |
|---|---|---|
| `/p5/sector <n>` | `-1` | 30-degree wedge 0..11; `-1` = all 12 |
| `/p5/layers <lo> <hi>` | `7 54` | TPC row range |
| `/p5/zrange <lo> <hi>` | `-105.5 105.5` | also `/p5/zrange full` (+-310) and `/p5/zrange tpc` |
| `/p5/pix/adcmin <v>` | `0` | drop pixels below this ADC |
| `/p5/pix/max <n>` | `250000` | cap on drawn pixels; stride-sampled above it, and it says so |

### Appearance — no re-read, instant

| command | default | meaning |
|---|---|---|
| `/p5/colour <mode>` | `level` | `level` / `class` / `track` (`layers` accepted as an alias for `level`) |
| `/p5/show <what> <bool>` | all true | `pix` / `clus` / `trk` / `g4` / `all` |
| `/p5/legend <bool>` | `true` | on-screen legend |
| `/p5/pix/adcmax <v>` | `1023` | ADC at the top of the colour ramp |
| `/p5/pix/size <px>` | `0` frame / `2.5` library | `0` = 1-px dots; `>0` = filled circles. The default differs by mode: a frame has ~200k pixels (dots are right), a library event ~5k G4 steps (dots vanish) |
| `/p5/pix/rowdr <bool>` | `true` | apply the v5.4 row-radius offsets to pixel r (false = GDML nominal) |
| `/p5/clussize <px>` | `2.2` | base screen size of a cluster mark (scaled up by cluster size) |

### Truth cluster chains ("tracks")

| command | default | meaning |
|---|---|---|
| `/p5/trk/minclus <n>` | `6` | minimum clusters for a chain |
| `/p5/trk/ptmin <GeV>` | `0.2` | minimum truth pT (above the 0.164 GeV looper wall) |
| `/p5/trk/primary <bool>` | `false` | keep only primary truth tracks |
| `/p5/trk/top <n>` | `3` | how many chains the `track` colour mode spotlights |
| `/p5/trk/select <id>` | `0` | spotlight one truth id; `0` clears it |
| `/p5/trk/list <n>` | — | print the N highest-pT chains (id, pT, primary, nclus, r range) |

### Geometry

| command | default | meaning |
|---|---|---|
| `/p5/geometry <preset>` | `tpc` | `tpc` (23/1020 volumes) / `tracker` (754) / `calo` (794) / `all` (1020) / `none` |
| `/p5/geom/show <substring>` | — | add every GDML volume whose name contains this |
| `/p5/geom/hide <substring>` | — | remove them |
| `/p5/geom/list [substring]` | — | print the names available to ask for |
| `/p5/geom/bright <f>` | `1.0` | geometry brightness. One map for both overlay modes — dim it here if a dense frame needs the geometry out of the way |
| `/p5/geom/tint <bool>` | `true` | per-system colouring; `false` = flat grey at the same luminance |
| `/p5/geom/style <s> [alpha]` | `auto` | `wireframe` / `surface [alpha]` / `auto`; alpha defaults 1.0 except 0.12 for surface |

A `/p5/geometry` preset resets everything, including earlier `geom/show` picks.

Each system carries its own tint and, if it has a hue, a legend row — uniform
grey is unreadable once more than the TPC is on screen, and brightness alone
does not fix that. The tints are dark, desaturated, and avoid the overlay's
blue/orange/aqua, so the data still wins.

| tier | system | volumes | colour |
|---|---|---|---|
| `tpc` | TPC gas | 1 | near-white (the subject) |
| `tpc` | TPC cage | 22 | mid grey |
| `tracker` | TPC endcap | 28 | dim grey |
| `tracker` | MVTX | 175 | magenta |
| `tracker` | INTT | 250 | violet |
| `tracker` | TPOT | 248 | yellow |
| `tracker` | beam pipe | 30 | grey |
| `calo` | CEMC | 3 | red |
| `calo` | HCal / cryostat | 6 | green |
| `calo` | EPD | 31 | dim grey |
| `all` | beam line | 226 | dimmest grey |

`/p5/geometry` prints this breakdown for whatever it just turned on.

### View and field

| command | default | meaning |
|---|---|---|
| `/p5/view <preset>` | `3d` | `3d` / `rphi` / `rz` / `iso` — viewpoint **and** framing |
| `/p5/zoom <radius_cm>` | per preset | frame this physical radius (TPC bounding radius ~131) |
| `/p5/field <file\|default\|off>` | off | load the field map so `/run/beamOn` tracks in field |

---

## 3. Standard Geant4 commands that matter here

Not part of `/p5/`, but this is a real Geant4 app, so they work — all verified in
this session:

| command | note |
|---|---|
| `/run/beamOn <n>` | tracks on the same geometry the batch g4hits came from; needs `/p5/field` for a realistic trajectory |
| `/gun/particle <name>`, `/gun/energy <E> <unit>`, `/gun/direction <x y z>` | only when no `--hepmc` was given |
| `/vis/viewer/set/style wireframe\|surface` | works now, but gives **opaque** chrome — prefer `/p5/geom/style surface`, which sets the alpha too |
| `/vis/viewer/set/viewpointThetaPhi <t> <p>` | `0 0` = r-phi, `90 0` = r-z |
| `/vis/viewer/zoomTo <f>`, `/vis/viewer/set/targetPoint <x y z> cm` | |
| `/vis/viewer/set/projection perspective <a> deg` | |
| `/vis/viewer/flush`, `/vis/viewer/rebuild` | flush writes the file in `TSG_OFFSCREEN` mode |
| `/vis/scene/endOfEventAction accumulate\|refresh` | for `/run/beamOn` trajectories |
| `/control/execute <macro>`, `/control/manual /p5/` | |
| `exit` | |

---

## 4. How dense should it look?

A **library event is one pp minimum-bias collision** and is sparse by physics,
not by setting. Measured over chunk 5 (1826 of 2000 generated events leave TPC
hits): median **2957** G4 steps/event, 10th pct 455, 90th pct 8922, max 24014 —
an order of magnitude of spread. If an event looks thin, it probably is one.

- `/p5/g4list 10` ranks the chunk's events and gives each as a multiple of the
  median; the menu prints the top 8 for you when the chunk has been indexed.
- Chunk 5 event 860 is 24014 hits = 8.1x median, and looks like a busy event.
- For a genuinely dense picture use **mode 1, a production frame**: ~200k
  pixels, ~20k clusters and ~250 truth chains of composed pile-up over the
  51 us readout. That is a different object, not a denser version of the same
  one — a frame is many collisions, a library event is exactly one.

## 4c. A frame is dense — what to turn down

The data palette is fine: at sector zoom the three layers separate cleanly
(blue pixel streaks, orange cluster dots, aqua chain lines). What defeats it is
volume — ~200k pixels over the whole TPC in a 1400 px window. No colour scheme
survives that, so look at less:

| knob | effect on frame 5 |
|---|---|
| `--sector 5` | 322,598 -> 16,788 px (one wedge of twelve) |
| `--adcmin 60` | 322,598 -> 77,582 px (drops the low-amplitude tail) |
| `--colour track` | everything dims except the leading truth chains |
| `--geometry tpc` | nothing but the gas and cage competing for attention |

The geometry colour map is the SAME in both overlay modes — a frame and a
library event show the same detector in the same colours, and the legend swatch
is the colour actually rendered (brightness and `/p5/geom/tint` included). If a
dense frame needs the detector out of the way, say so explicitly with
`/p5/geom/bright 0.5` or a lighter `/p5/geometry` preset.

Reads well: `./p5_display.sh --frame 5 --adcmin 60 --geometry tracker`

## 4b. Recipes

```bash
# the frame that is actually on screen, all 12 sectors
./p5_display.sh --frame 12

# ML label view of one wedge, tracker geometry
./p5_display.sh --frame 12 --sector 5 --colour class --geometry tracker

# the whole detector, legible: everything visible needs brighter chrome
./p5_display.sh --g4chunk 8 --g4event 2 --geometry all --geombright 2.2

# translucent shells, beam-view, straight from the shell
./p5_display.sh --frame 12 --geomstyle surface 0.10 --view rphi

# spotlight one truth chain, no pixels
./p5_display.sh --frame 12 --cmd "/p5/trk/select 5803" --cmd "/p5/show pix false"

# headless still, no GUI at all
./p5_display.sh --batch --vis TSG_OFFSCREEN --frame 12 --colour class \
    --view rphi --cmd "/vis/viewer/flush"

# the raw g4 steps of one library collision
./p5_display.sh --g4chunk 0 --g4event 7 --view rphi

# a BUSY library collision (find one with /p5/g4list, then load it)
./p5_display.sh --g4chunk 5 --g4event 860 --geometry tracker
```
