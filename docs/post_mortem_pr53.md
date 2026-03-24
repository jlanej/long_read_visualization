# Post Mortem — PR #53: Add softclips/mismatches toggles, lock hap panel navigation

**PR:** https://github.com/jlanej/long_read_visualization/pull/53  
**Merged:** 2026-03-24  
**Reverted by:** this PR  

---

## Summary

PR #53 introduced two features: display-toggle buttons for soft clips and
mismatches (alongside a hardening of the existing indel toggle), and a
`_lockHapPanel()` function that blocked mouse/touch/wheel navigation on the
Haplotype 1 and Haplotype 2 IGV panels.  The intended goal was to prevent
accidental drag-to-pan or scroll-zoom that would de-sync or unload reads in the
hap panels.  Both features shipped with tests that passed; nonetheless the PR was
reverted because the implementation contained fundamental UI defects that the
tests were not equipped to catch.

---

## What Was Changed

| File | Change |
|------|--------|
| `server/static/index.html` | Added `toggleSoftClips` and `toggleMismatches` toolbar buttons; added `showSoftClips`/`showMismatches` state variables; replaced hard-coded `showSoftClips: false` in track configs with the state variable; added `showMismatches` to track configs; added `updateAllTrackSoftClipDisplay()`, `updateAllTrackMismatchDisplay()`; added `_lockHapPanel()` and called it for both hap panels after browser creation |
| `tests/test_server.py` | Added `TestFrontendDisplayToggles` (14 tests) and `TestFrontendHapPanelLockdown` (12 tests) |
| `docs/visualization_guide.md` | Added soft clips / mismatches toggle entries; added hap panel navigation-lock note |
| `scripts/generate_docs.py` | Updated doc-generation strings to match the above |

---

## Root Causes of Failure

### 1. `_lockHapPanel` — Wheel handler incorrectly assumed a CSS-scroll container

**The code:**
```js
el.addEventListener("wheel", (e) => {
    let scroller = e.target;
    while (scroller && scroller !== el) {
        if (scroller.scrollHeight > scroller.clientHeight) {
            scroller.scrollTop += e.deltaY;   // forward vertical scroll
            break;
        }
        scroller = scroller.parentElement;
    }
    e.preventDefault();   // block IGV.js zoom
    e.stopPropagation();
}, { capture: true, passive: false });
```

**Why it fails:** IGV.js renders alignment reads on a `<canvas>` element using
an *internal virtual-scroll* model.  It tracks its own `scrollTop` state and
redraws via `requestAnimationFrame`; the DOM containers that wrap the canvas use
`overflow: hidden`, not `overflow: scroll` or `overflow: auto`.  As a result
`scrollHeight` equals `clientHeight` for every element in the chain and the
`while` loop terminates without ever calling `scroller.scrollTop += e.deltaY`.
The forwarded vertical scroll is silently discarded and **vertical scrolling
through stacked reads stops working entirely** — both with a mouse wheel and a
trackpad.

### 2. `_lockHapPanel` — `touchmove` was unconditionally blocked

**The code:**
```js
el.addEventListener("touchmove", (e) => {
    e.preventDefault();
    e.stopPropagation();
}, { capture: true, passive: false });
```

**Why it fails:** This intercepts *every* touch-move event in the panel,
including two-finger vertical scroll.  On touch screens and tablets users
cannot scroll down through a tall pile-up of reads at all.  There is no
equivalent "forward vertical touch scroll" logic, unlike the (already broken)
wheel handler.

### 3. `_lockHapPanel` — `pointermove` suppression introduced drag latency and broke hover

**The code:**
```js
el.addEventListener("pointermove", (e) => {
    if (dragOriginX !== null &&
        Math.abs(e.clientX - dragOriginX) > DRAG_PX) {
        e.stopPropagation();
    }
}, { capture: true });
```

**Why it fails:** IGV.js listens to `pointermove` in its own capture-phase
handlers for read hovering, tooltip display, and cursor-shape updates.
`stopPropagation()` at the capture phase terminates the event chain *before*
any IGV.js handler sees the event, so tooltips and hover highlights stop
working once the user has moved the pointer more than 3 px horizontally.
The 3 px threshold is also too tight — normal read-selection clicks involving
a tiny tremor can exceed it.

### 4. Display toggle update functions used an ineffective IGV.js API

**The code:**
```js
function updateAllTrackSoftClipDisplay() {
    for (const browser of [refBrowser, hap1Browser, hap2Browser]) {
        if (!browser) continue;
        for (const tv of browser.trackViews) {
            const t = tv.track;
            if (t && t.type === "alignment") {
                t.showSoftClips = showSoftClips;
                if (t.config) t.config.showSoftClips = showSoftClips;
            }
        }
        browser.updateViews();
    }
}
```

**Why it fails:** Mutating `track.showSoftClips` and `track.config.showSoftClips`
directly bypasses IGV.js's internal state management.  The `browser.updateViews()`
call redraws the current viewport but IGV.js re-reads alignment rendering
options at the time it queries the reference genome for the current locus —
the mutated properties may be overwritten by a pending async fetch callback.
In practice the toggle buttons appear to work in isolation but revert silently
when the browser pans or zooms, creating a confusing, unreliable experience.

### 5. Tests were limited to static HTML string matching

All 26 new tests in `TestFrontendDisplayToggles` and `TestFrontendHapPanelLockdown`
asserted that specific strings (function names, variable names, event listener
calls) are present in `index.html`.  They confirmed code *exists* but did not
validate runtime behaviour.  A Playwright or Selenium integration test running
against a real IGV.js instance would have caught defects 1–4 immediately.

---

## Evidence That This PR Reverts All PR #53 Changes

### Diff summary (revert branch vs. main / post-PR-53)

| File | Lines removed | Lines added |
|------|--------------|-------------|
| `server/static/index.html` | 278 | 2 |
| `tests/test_server.py` | 143 | 0 |
| `docs/visualization_guide.md` | 6 | 0 |
| `scripts/generate_docs.py` | 9 | 0 |

These numbers are the **exact inverse** of the PR #53 diff
(`+278 -2`, `+143 -0`, `+6 -0`, `+9 -0` respectively), confirming a
complete revert.

### Specific patterns removed

The following identifiers introduced by PR #53 are **absent** from every
file after this revert:

- `toggleSoftClips`, `showSoftClips`, `updateAllTrackSoftClipDisplay`
- `toggleMismatches`, `showMismatches`, `updateAllTrackMismatchDisplay`
- `_lockHapPanel`, `DRAG_PX`, `dragOriginX`, `scroller.scrollTop += e.deltaY`
- `TestFrontendDisplayToggles`, `TestFrontendHapPanelLockdown`
- soft clips / mismatches / navigation-lock text in docs and doc generator

Verified by:
```
grep -rn "toggleSoftClips\|showSoftClips\|_lockHapPanel\|DRAG_PX\|TestFrontendDisplayToggles" \
     server/static/index.html tests/test_server.py docs/ scripts/generate_docs.py
# → no matches
```

### Test suite

After this revert, **368 tests pass, 33 skipped** (identical pass rate to the
baseline before PR #53 was opened), with zero failures or errors.

---

## Lessons Learned

1. **Don't intercept low-level pointer/touch events inside an embedded
   third-party rendering engine.**  IGV.js owns its own event model; fighting
   it with capture-phase `stopPropagation` causes unpredictable regressions.

2. **Never assume CSS `overflow` properties for a library's internal scroll.**
   Verify the actual DOM structure with DevTools before writing a scroll-
   forwarding shim.

3. **Static HTML string tests do not exercise runtime behaviour.**  Any future
   UI change that affects user interaction (dragging, scrolling, toggling)
   should be validated by an end-to-end browser test (e.g. Playwright).

4. **Prefer IGV.js public API for property changes** (`browser.loadTrack`,
   `track.setConfiguration`, etc.) over direct property mutation followed by
   `updateViews()`, which can race against pending async renders.
