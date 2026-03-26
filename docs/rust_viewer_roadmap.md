# Rust Viewer Roadmap: Feature Parity and Beyond

> **Meta Issue**: Bring the Rust native viewer (`viewer/`) to full feature parity with
> the Python + igv.js implementation (`server/`), then exceed it with native-app
> advantages (performance, offline use, richer interactivity).
> Ensure that all features are well tested via automated CI, and ensure that accuracy and interpretation of the visualization is considered with each implementation. 
> Please ensure that our methods you come across, create, or edit are reasonably efficient and not something like O(n squared) for very large compute tasks that can annoying block interactivity. Keep performance in mind for all implementations 

## Background

The Rust viewer was ported from the Python/igv.js web UI in issues #56–#62. Those
issues established the foundational architecture:

- Three-panel layout (Reference, Haplotype 1, Haplotype 2)
- BAM/CRAM + FASTA file I/O via `noodles`
- CIGAR-aware coordinate translation
- Basic pileup rendering with HP tag coloring
- Region navigation from manifests and VCFs
- K-mer dot plots

However, a careful comparison of the two implementations reveals **substantial
missing functionality** and **performance gaps** in the Rust viewer. This roadmap
enumerates every gap, organizes them into actionable sub-issues, and prioritizes
them so that the viewer becomes a production-quality genomics review tool.

---

## Gap Summary

| Category | Python/igv.js | Rust Viewer | Status |
|----------|--------------|-------------|--------|
| Alignment tracks per panel | 3 (reads + 2 assembly cross-alignments) | 1 (reads only) | ❌ Missing |
| Mismatch/SNV coloring | Colored bases at mismatch positions | Not rendered | ❌ Missing |
| Soft-clip display | Toggle-able soft-clip bases | Parsed but not rendered | ❌ Missing |
| Indel visual symbols | igv.js native indel markers | Overlaid rectangles only | ⚠️ Partial |
| Read sorting by HP tag | Reads grouped/sorted by haplotype | Sorted by start position only | ❌ Missing |
| Hover tooltips | igv.js built-in read popups | `read_name` stored but no UI | ❌ Missing |
| Coverage visualization | Implicit via read density | No coverage track | ❌ Missing |
| Async data loading | Threaded HTTP server + igv.js cache | Synchronous, blocks UI thread | ❌ Missing |
| Data caching | LRU cache (256) + igv.js internal | No caching, full reload each nav | ❌ Missing |
| Dynamic coord sync | Debounced re-translation on every pan/zoom | Static: region-load only | ⚠️ Partial |
| Sample selector | Multi-sample dropdown | First sample only | ❌ Missing |
| Base-level sequence view | igv.js at high zoom shows letters | First 40 bp text summary | ❌ Missing |
| Mouse pan/zoom | igv.js native drag + scroll | Buttons only | ❌ Missing |
| Coordinate ruler | igv.js genomic axis | None | ❌ Missing |
| Compare Reads modal | Cross-panel read ID matching | None | ❌ Missing |
| Debounce/throttle | 150–350 ms debounced sync | None | ❌ Missing |
| Default indel threshold | 50 bp (long-read optimized) | 3 bp (short-read default) | ⚠️ Mismatch |

---

## Sub-Issues

### Phase 1 — Critical Feature Parity (Visualization Accuracy)

These issues address the most impactful visual and data-accuracy gaps.
Without them, the Rust viewer cannot be used for reliable genomic review.

#### Sub-Issue 1: Load and Display Cross-Alignment (Assembly) Tracks

**Priority:** P0 — Blocks meaningful use  
**Effort:** Large  
**Files:** `viewer/src/app.rs`, `viewer/src/genome/bam.rs`, `viewer/src/genome/pileup.rs`

**Problem:**  
The Python/igv.js version shows up to **9 alignment tracks** across 3 panels:

| Panel | Reads Track | Assembly Track 1 | Assembly Track 2 |
|-------|-------------|-------------------|-------------------|
| Reference | Reads → Ref | Hap1 → Ref | Hap2 → Ref |
| Haplotype 1 | Reads → Hap1 | Ref → Hap1 | Hap2 → Hap1 |
| Haplotype 2 | Reads → Hap2 | Ref → Hap2 | Hap1 → Hap2 |

The Rust viewer loads **only 3 tracks** (reads per panel). Assembly cross-alignments
(hap-to-ref, ref-to-hap, hap-to-hap) are completely absent.

**Requirements:**
- Extend `DataPaths` to include paths for all 6 cross-alignment BAM files
  (hap1_to_ref, hap2_to_ref, ref_to_hap1, ref_to_hap2, hap1_to_hap2, hap2_to_hap1)
- Extend `from_tsv()` file discovery to locate these BAMs from preprocessing output
- Add per-panel secondary track rendering (stacked below reads track)
- Color assembly tracks distinctly (matching Python: blue for ref, green for hap1, orange for hap2)
- Assembly tracks default to squished display mode

**Acceptance Criteria:**
- All 9 tracks visible when BAMs are available
- Assembly tracks use distinct color per source assembly
- Graceful fallback when cross-alignment BAMs are absent
- Tests for multi-track loading and rendering

---

#### Sub-Issue 2: Render Mismatches (SNVs) as Colored Bases on Reads

**Priority:** P0 — Critical for variant review  
**Effort:** Medium  
**Files:** `viewer/src/genome/bam.rs`, `viewer/src/genome/pileup.rs`

**Problem:**  
igv.js displays individual mismatched bases as colored letters/blocks on reads,
making SNVs immediately visible. The Rust viewer renders reads as solid colored
rectangles with no base-level detail.

The `AlignedRead` struct already has a `mismatches` field (genome/bam.rs) but it
is never populated or rendered.

**Requirements:**
- During BAM/CRAM query, walk the CIGAR string + reference sequence to identify mismatches
- Store mismatch positions and bases in `AlignedRead.mismatches`
- In `layout_read_rects()`, overlay colored markers at mismatch positions
- Use standard genomics base colors (A=green, T=red, C=blue, G=orange/yellow)
- Only render mismatches when zoom level shows ≥1 px per base
- Add toolbar toggle ("Show Mismatches") matching Python behavior

**Acceptance Criteria:**
- Mismatches visible at sufficient zoom levels
- Base colors follow genomics convention
- Toggle works to show/hide mismatches
- Tests for mismatch extraction from CIGAR + reference

---

#### Sub-Issue 3: Display Soft-Clipped Bases

**Priority:** P1 — Important for SV breakpoint review  
**Effort:** Small  
**Files:** `viewer/src/genome/bam.rs`, `viewer/src/genome/pileup.rs`

**Problem:**  
Soft clips are parsed from BAM CIGAR strings but never rendered. For SV review,
soft-clipped bases are critical indicators of breakpoint positions and novel
insertions.

**Requirements:**
- Extend read rendering to show soft-clipped regions as semi-transparent extensions
- Use distinct color/opacity for clipped bases (e.g., lighter shade of read color)
- Add toolbar toggle ("Show Soft Clips") matching Python behavior (default: off)
- Soft clips should extend read boundaries appropriately

**Acceptance Criteria:**
- Soft clips visible when toggle is enabled
- Default off (matching Python)
- Clipped region visually distinct from aligned region
- Tests for soft-clip rendering

---

#### Sub-Issue 4: Sort and Group Reads by Haplotype (HP Tag)

**Priority:** P1 — Important for phasing review  
**Effort:** Small  
**Files:** `viewer/src/genome/pileup.rs`

**Problem:**  
Python/igv.js sorts reads by HP tag (hap1 reads grouped together, hap2 together,
unphased at bottom). The Rust viewer only sorts by start position, causing hap1
and hap2 reads to interleave, making phasing patterns hard to see.

**Requirements:**
- Add HP-tag-based sorting as a packing mode in `pack_reads()`
- Group order: HP=1 (top), HP=2 (middle), unphased (bottom)
- Within each group, maintain start-position sorting
- Add toolbar toggle for "Sort by HP" vs "Sort by Position"

**Acceptance Criteria:**
- When enabled, reads visually separate into haplotype groups
- Clear visual boundary between HP groups
- Tests for HP-sorted packing with mixed haplotype reads

---

#### Sub-Issue 5: Correct Default Indel Threshold for Long Reads

**Priority:** P1 — Quick fix, high impact  
**Effort:** Trivial  
**Files:** `viewer/src/genome/pileup.rs`

**Problem:**  
The Rust viewer defaults to hiding indels ≤3 bp (`DEFAULT_INDEL_THRESHOLD = 3`),
which is a short-read setting. The Python/igv.js version uses 50 bp as the
default threshold, appropriate for long-read data where small indels are noisy.

**Requirements:**
- Change `DEFAULT_INDEL_THRESHOLD` from 3 to 50
- Ensure the toolbar threshold spinner reflects the new default
- Add a "3rd Gen View" or "Long Read Mode" preset that sets optimal defaults

**Acceptance Criteria:**
- Default threshold is 50 bp on fresh launch
- Users can still adjust threshold manually
- Test updated default value

---

### Phase 2 — Performance and Responsiveness

These issues address the performance gaps that make the viewer
sluggish or unusable for production-size datasets.

#### Sub-Issue 6: Async/Threaded Data Loading (Non-Blocking UI)

**Priority:** P0 — Blocks production use  
**Effort:** Large  
**Files:** `viewer/src/app.rs`, new module `viewer/src/data_loader.rs`

**Problem:**  
All BAM/CRAM queries and FASTA lookups in `load_region_data()` execute
synchronously on the main UI thread. For large BAM files or regions deep into
chromosomes (issue #33), this freezes the UI for seconds to minutes.

**Requirements:**
- Create a background data-loading thread/channel architecture
- `load_region_data()` sends a request to background thread and returns immediately
- Background thread performs BAM queries, FASTA lookups, pileup packing
- Results are sent back to UI thread via channel (e.g., `std::sync::mpsc`)
- UI shows a loading spinner/indicator while data loads
- Support cancellation: if user navigates before load completes, cancel stale request

**Implementation Guidance:**
- Use `std::thread::spawn` + `mpsc::channel` pattern
- egui's `ctx.request_repaint()` to wake UI when data arrives
- Consider `Arc<AtomicBool>` cancellation token pattern

**Acceptance Criteria:**
- UI remains responsive during BAM loading
- Loading indicator visible during queries
- Stale loads are cancelled on new navigation
- Tests for concurrent load/cancel behavior

---

#### Sub-Issue 7: Region Data Caching (LRU Cache)

**Priority:** P1 — Important for navigation UX  
**Effort:** Medium  
**Files:** `viewer/src/app.rs`

**Problem:**  
Every region navigation triggers a full reload of all BAM/FASTA data. When
reviewing SVs, users frequently navigate back and forth between regions. The
Python version uses LRU caching (256 entries) to avoid redundant loads.

**Requirements:**
- Implement LRU cache for loaded panel data (keyed by region string)
- Cache should store packed pileup rows + FASTA sequences
- Configurable cache size (default: 32 regions)
- Cache invalidation on file path changes

**Acceptance Criteria:**
- Revisiting a previously loaded region is near-instant
- Cache eviction works correctly (LRU order)
- Memory usage is bounded
- Tests for cache hit/miss/eviction behavior

---

#### Sub-Issue 8: Dynamic Coordinate Translation During Pan/Zoom

**Priority:** P1 — Important for interactive exploration  
**Effort:** Medium  
**Files:** `viewer/src/panel_sync.rs`, `viewer/src/app.rs`

**Problem:**  
The Python version re-translates coordinates on every pan/zoom event (with 150–350
ms debounce). The Rust viewer only computes haplotype regions at region-load time.
When the user manually pans or zooms a panel, the other panels drift out of sync.

**Requirements:**
- On pan/zoom events, re-run coordinate translation for the new view window
- Implement debounce (150 ms for proportional, 350 ms for index-based mapping)
- Use the async loading infrastructure (Sub-Issue 6) for non-blocking translation
- Update panel views with translated coordinates

**Acceptance Criteria:**
- Panels remain synchronized during interactive pan/zoom
- Debounce prevents excessive re-computation
- Works with both proportional and index-based mapping
- Tests for sync accuracy after pan/zoom sequences

---

### Phase 3 — Modern UI and Interactivity

These issues bring the viewer to a modern, professional standard
and add capabilities that exceed the web-based version.

#### Sub-Issue 9: Mouse-Driven Pan and Zoom

**Priority:** P0 — Expected in any genome browser  
**Effort:** Medium  
**Files:** `viewer/src/app.rs`, `viewer/src/panel_sync.rs`

**Problem:**  
The Rust viewer only supports pan/zoom via toolbar buttons and keyboard shortcuts.
Every modern genome browser (igv.js, IGV desktop, JBrowse) supports click-drag
panning and scroll-wheel zooming.

**Requirements:**
- Click-drag on a panel pans the view (horizontal drag → genomic coordinate shift)
- Scroll wheel zooms in/out centered on cursor position
- Respect sync mode: if sync enabled, propagate to other panels
- Use `PanelSyncManager::view_mut()` (currently dead code in issue #72) for direct panel manipulation
- Smooth animation/interpolation for zoom transitions

**Acceptance Criteria:**
- Drag-to-pan works on all three panels
- Scroll-to-zoom works with cursor-centered zoom
- Sync propagation works correctly
- Tests for mouse event → coordinate translation

---

#### Sub-Issue 10: Hover Tooltips with Read Metadata

**Priority:** P1 — Important for variant review  
**Effort:** Medium  
**Files:** `viewer/src/genome/pileup.rs`, `viewer/src/app.rs`

**Problem:**  
igv.js provides rich tooltips on read hover/click (read name, MAPQ, strand, flags,
HP tag, base quality). The Rust viewer stores `read_name` in `ReadRect` but never
displays it. Issue #72 lists `ReadRect::read_name` as dead code needing tooltip UI.

**Requirements:**
- Detect mouse hover over read rectangles in pileup
- Display egui tooltip with:
  - Read name
  - Mapping quality
  - HP tag / haplotype assignment
  - Strand (forward/reverse)
  - Start–end coordinates
  - CIGAR summary (alignment length, # indels)
- Tooltip should not obscure the view significantly

**Acceptance Criteria:**
- Hovering over a read shows metadata tooltip
- Tooltip content is accurate
- Performance: tooltip lookup is O(1) or O(log n), not O(n) scan
- Tests for tooltip content generation

---

#### Sub-Issue 11: Genomic Coordinate Ruler / Scale Bar

**Priority:** P1 — Important for orientation  
**Effort:** Medium  
**Files:** `viewer/src/app.rs`, new module `viewer/src/ruler.rs`

**Problem:**  
igv.js shows a genomic coordinate axis with tick marks and position labels. The
Rust viewer has no coordinate ruler, making it difficult to determine exact
positions. Issue #72 lists `FastaSequence::{start, end}` as dead code needing
a coordinate ruler.

**Requirements:**
- Render a coordinate ruler at the top of each panel
- Show chromosome name + position range
- Tick marks at readable intervals (auto-scaled: 1bp, 10bp, 100bp, 1kb, etc.)
- Position labels at major tick marks
- Ruler updates with pan/zoom

**Acceptance Criteria:**
- Ruler shows accurate coordinates at all zoom levels
- Tick spacing auto-adjusts for readability
- Ruler is visually clean and non-intrusive
- Tests for tick interval calculation

---

#### Sub-Issue 12: Base-Level Sequence Display at High Zoom

**Priority:** P1 — Important for detailed variant review  
**Effort:** Medium  
**Files:** `viewer/src/app.rs`, `viewer/src/genome/pileup.rs`

**Problem:**  
Currently the Rust viewer shows only a text summary of the first 40 bp of FASTA
sequence. igv.js renders individual nucleotide letters when zoomed in far enough.
For reviewing individual variants, seeing the actual sequence is essential.

**Requirements:**
- When the zoom level reaches ≥ ~5 px per base, render individual nucleotide letters
- Use standard base colors (A=green, T=red, C=blue, G=orange)
- Show reference sequence as a track above reads in each panel
- Show read-level bases when at very high zoom (≥ ~10 px per base)
- Align read bases to reference for mismatch visualization

**Acceptance Criteria:**
- Sequence letters appear at high zoom
- Colors follow standard genomics convention
- Performance remains smooth during zoom transitions
- Tests for sequence rendering thresholds

---

#### Sub-Issue 13: Coverage Track / Depth Histogram

**Priority:** P2 — Nice to have for overview  
**Effort:** Medium  
**Files:** `viewer/src/genome/pileup.rs`, `viewer/src/app.rs`

**Problem:**  
igv.js shows read density as coverage. The Rust viewer has no coverage
visualization, making it hard to spot dropout regions or coverage imbalances.

**Requirements:**
- Compute per-base (or per-bin) coverage from loaded reads
- Render as a histogram above the pileup in each panel
- Height proportional to coverage depth
- Optional HP-stratified coverage (stacked: hap1 green + hap2 orange)
- Configurable height (collapsible)

**Acceptance Criteria:**
- Coverage histogram visible above reads
- Scales correctly with zoom/pan
- HP-stratified view works
- Tests for coverage computation

---

#### Sub-Issue 14: Multi-Sample Support with Sample Selector

**Priority:** P2 — Important for production workflows  
**Effort:** Medium  
**Files:** `viewer/src/app.rs`, `viewer/src/main.rs`

**Problem:**  
The Python version supports a sample dropdown to switch between multiple samples
in a TSV config. The Rust viewer loads only the first sample (line 91 in app.rs).

**Requirements:**
- Parse all samples from TSV config (not just first)
- Add sample selector dropdown in toolbar
- On sample switch: update all file paths and reload data
- Maintain region position across sample switches when possible

**Acceptance Criteria:**
- All samples from TSV are selectable
- Switching samples correctly updates all panels
- Region position preserved if same manifest loaded
- Tests for multi-sample config parsing

---

#### Sub-Issue 15: Compare Reads Across Panels

**Priority:** P2 — Useful for phasing validation  
**Effort:** Medium  
**Files:** `viewer/src/app.rs`

**Problem:**  
The Python version has a "Compare Reads" modal that finds reads present in
multiple panels (i.e., reads mapped to both reference and haplotype). This helps
validate phasing accuracy. The Rust viewer has no equivalent.

**Requirements:**
- Collect read names from all loaded panels
- Compute intersection/union statistics (reads in all 3, in 2, in only 1)
- Display summary in a modal/side panel
- Optionally highlight shared reads in pileup

**Acceptance Criteria:**
- Read comparison statistics are accurate
- Modal displays overlap counts
- Shared reads can be visually identified
- Tests for set intersection computation

---

### Phase 4 — Polish and Production Readiness

#### Sub-Issue 16: SV Annotation Overlay from Coordinate Translation

**Priority:** P2  
**Effort:** Medium  
**Files:** `viewer/src/app.rs`, `viewer/src/genome/coordinate_mapper.rs`

**Problem:**  
Issue #72 lists `QueryResult` fields (ref_chrom, strand, mapq, gap sizes) as dead
code needing an SV annotation overlay. The coordinate mapper classifies SV gaps
(deletion, insertion, inversion, translocation, complex) but this information is
never displayed in the UI.

**Requirements:**
- When loading a region, run coordinate translation and collect SV gap events
- Render SV annotations as labeled markers between panels or in panel margins
- Show gap type, size, and strand information
- Color-code by SV type (deletions=red, insertions=blue, inversions=purple, etc.)

**Acceptance Criteria:**
- SV gap annotations visible between panels
- Gap types correctly labeled
- Annotations update on region navigation
- Tests for SV annotation rendering

---

#### Sub-Issue 17: Variant Metadata Display in UI

**Priority:** P2  
**Effort:** Small  
**Files:** `viewer/src/app.rs`, `viewer/src/region.rs`

**Problem:**  
Issue #72 lists manifest serde fields (genotype, fasta_region, description) as
dead code. When navigating regions from a manifest, the variant metadata
(genotype, SV type, size) should be displayed prominently.

**Requirements:**
- Show variant metadata (genotype, size, type) in toolbar or panel header
- Display region description if available
- Show coordinate mapping information (ref → hap translation result)

**Acceptance Criteria:**
- Variant metadata visible during region review
- All manifest fields rendered
- Tests for metadata display

---

#### Sub-Issue 18: Export / Screenshot Functionality

**Priority:** P3  
**Effort:** Small  
**Files:** `viewer/src/app.rs`

**Problem:**  
No way to save the current view for reports or publications.

**Requirements:**
- Add "Export PNG" button to save current viewport
- Include all panels, ruler, and annotations in export
- Use egui's screenshot or custom rendering path

**Acceptance Criteria:**
- Exported PNG is pixel-accurate
- All visual elements included
- File picker for save location

---

#### Sub-Issue 19: Accessibility — Colorblind-Friendly Palette

**Priority:** P2  
**Effort:** Small  
**Files:** `viewer/src/genome/pileup.rs`, `viewer/src/app.rs`

**Problem:**  
Current HP tag colors (green/orange) are not optimal for colorblind users
(red-green colorblindness affects ~8% of males). No alternative palette is offered.

**Requirements:**
- Add colorblind-safe palette option (e.g., blue/orange from ColorBrewer)
- Settable via toolbar or preferences
- Apply to all colored elements (reads, assembly tracks, dot plot)

**Acceptance Criteria:**
- Colorblind palette clearly distinguishes haplotypes
- Toggle is easily accessible
- Tests for palette color values

---

#### Sub-Issue 20: Performance Profiling and Memory Optimization (Issue #33)

**Priority:** P1 — Blocks real-world use  
**Effort:** Large  
**Files:** `viewer/src/genome/bam.rs`, `viewer/src/genome/coordinate_mapper.rs`

**Problem:**  
Issue #33 reports memory overload and performance issues for regions deep into
chromosomes. Regions near chromosome starts load quickly; those farther in take
a long time or fail with OOM.

**Requirements:**
- Profile BAM index seeking for mid-chromosome regions
- Profile coordinate mapper index loading and query for large chromosomes
- Identify and fix O(n) scans or excessive memory allocation
- Add benchmarks for region loading at various chromosome positions
- Consider streaming/lazy loading for large coordinate indices
- Implement the memory-bounded data loading from Sub-Issue 7

**Acceptance Criteria:**
- Regions at all chromosome positions load within 5 seconds
- Memory usage stays under 2 GB for typical SV review
- Benchmark suite validates performance
- Issue #33 resolved

---

## Implementation Order

```
Phase 1 (Feature Parity — Weeks 1–4):
  ├── Sub-Issue 5:  Fix indel threshold default .............. [Trivial, Day 1]
  ├── Sub-Issue 4:  HP tag sorting ........................... [Small, Week 1]
  ├── Sub-Issue 2:  Mismatch rendering ....................... [Medium, Week 1–2]
  ├── Sub-Issue 3:  Soft-clip display ........................ [Small, Week 2]
  └── Sub-Issue 1:  Cross-alignment tracks ................... [Large, Week 2–4]

Phase 2 (Performance — Weeks 3–6):
  ├── Sub-Issue 6:  Async data loading ....................... [Large, Week 3–5]
  ├── Sub-Issue 20: Performance profiling (#33) .............. [Large, Week 4–6]
  ├── Sub-Issue 7:  LRU caching .............................. [Medium, Week 5]
  └── Sub-Issue 8:  Dynamic coordinate sync .................. [Medium, Week 5–6]

Phase 3 (Modern UI — Weeks 5–8):
  ├── Sub-Issue 9:  Mouse pan/zoom ........................... [Medium, Week 5]
  ├── Sub-Issue 10: Hover tooltips ........................... [Medium, Week 6]
  ├── Sub-Issue 11: Coordinate ruler ......................... [Medium, Week 6]
  ├── Sub-Issue 12: Base-level sequence ...................... [Medium, Week 7]
  └── Sub-Issue 13: Coverage track ........................... [Medium, Week 7–8]

Phase 4 (Polish — Weeks 7–10):
  ├── Sub-Issue 14: Multi-sample support ..................... [Medium, Week 7]
  ├── Sub-Issue 15: Compare reads modal ...................... [Medium, Week 8]
  ├── Sub-Issue 16: SV annotation overlay .................... [Medium, Week 8]
  ├── Sub-Issue 17: Variant metadata display ................. [Small, Week 9]
  ├── Sub-Issue 18: Export/screenshot ........................ [Small, Week 9]
  └── Sub-Issue 19: Colorblind palette ....................... [Small, Week 9]
```

## Relationship to Existing Issues

| Existing Issue | Related Sub-Issues | Notes |
|---------------|-------------------|-------|
| #72 (Dead code items) | 9, 10, 11, 12, 16, 17 | Dead code items become live as features are implemented |
| #33 (Memory/perf) | 6, 7, 20 | Directly addressed by async loading, caching, profiling |
| #14 (Old TODO) | 1, 6, 5 | CRAM, read overlap tracking, sync, indel defaults |
| #38–#44 (Graph genomes) | Future | Graph features build on top of the foundation established here |

## Success Criteria

The Rust viewer is considered at feature parity when:

1. **All 9 alignment tracks** load and render correctly per panel
2. **Mismatches, indels, and soft clips** are visually rendered on reads
3. **Reads are HP-sorted** with clear haplotype grouping
4. **UI remains responsive** during all data loading operations
5. **Mouse pan/zoom** works natively in all panels
6. **Hover tooltips** show read metadata on interaction
7. **Coordinate ruler** shows genomic positions
8. **Base-level sequence** renders at high zoom
9. **Multiple samples** can be loaded and switched
10. **Performance**: any region loads within 5 seconds, memory under 2 GB

The viewer exceeds the web version when:

1. **Offline operation** with no server dependency
2. **Sub-frame rendering latency** for smooth interaction
3. **Native OS integration** (file dialogs, keyboard shortcuts, clipboard)
4. **Dot plot** computation is instant (no server round-trip)
5. **Cross-platform binary** distribution via CI/CD
