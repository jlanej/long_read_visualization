use std::collections::HashMap;
use std::collections::VecDeque;

use crate::genome::FastaSequence;
use crate::genome::pileup::PileupRow;

// ---------------------------------------------------------------------------
// Cached panel data
// ---------------------------------------------------------------------------

/// Assembly cross-alignment track data for caching.
#[derive(Debug, Clone)]
pub struct CachedAssemblyTrack {
    pub label: String,
    pub color: [u8; 3],
    pub rows: Vec<PileupRow>,
}

/// All data for one panel, stored in the cache.
#[derive(Debug, Clone, Default)]
pub struct CachedPanelData {
    pub rows: Vec<PileupRow>,
    pub sequence: Option<FastaSequence>,
    pub assembly_tracks: Vec<CachedAssemblyTrack>,
}

/// A complete set of panel data for one region navigation event.
#[derive(Debug, Clone)]
pub struct CachedRegion {
    pub ref_data: CachedPanelData,
    pub hap1_data: CachedPanelData,
    pub hap2_data: CachedPanelData,
    pub status_message: String,
}

// ---------------------------------------------------------------------------
// LRU Cache
// ---------------------------------------------------------------------------

/// Default maximum number of cached regions.
pub const DEFAULT_CACHE_CAPACITY: usize = 32;

/// LRU cache for loaded region data, keyed by a region string.
///
/// Uses a `HashMap` for O(1) lookup and a `VecDeque<String>` for LRU ordering.
/// The most recently used key is at the back of the deque.
pub struct RegionCache {
    map: HashMap<String, CachedRegion>,
    order: VecDeque<String>,
    capacity: usize,
}

impl RegionCache {
    /// Create a new cache with the given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            map: HashMap::new(),
            order: VecDeque::with_capacity(capacity),
            capacity,
        }
    }

    /// Look up a region in the cache.  If found, promote it to most-recently
    /// used and return a reference.
    pub fn get(&mut self, key: &str) -> Option<&CachedRegion> {
        if self.map.contains_key(key) {
            self.promote(key);
            self.map.get(key)
        } else {
            None
        }
    }

    /// Insert (or update) a region in the cache.  If the cache is at capacity,
    /// the least-recently used entry is evicted first.
    pub fn insert(&mut self, key: String, value: CachedRegion) {
        if self.map.contains_key(&key) {
            self.promote(&key);
            self.map.insert(key, value);
            return;
        }

        if self.order.len() >= self.capacity {
            // Evict the least-recently used (front of the deque) — O(1).
            if let Some(evicted) = self.order.pop_front() {
                self.map.remove(&evicted);
            }
        }
        self.order.push_back(key.clone());
        self.map.insert(key, value);
    }

    /// Clear the entire cache (e.g. on file path change).
    pub fn clear(&mut self) {
        self.map.clear();
        self.order.clear();
    }

    /// Number of entries currently in the cache.
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// Whether the cache is empty.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// Maximum capacity.
    #[allow(dead_code)]
    pub fn capacity(&self) -> usize {
        self.capacity
    }

    // -- internal -----------------------------------------------------------

    /// Move `key` to the most-recently-used position (back of `order`).
    fn promote(&mut self, key: &str) {
        if let Some(pos) = self.order.iter().position(|k| k == key) {
            self.order.remove(pos);
            self.order.push_back(key.to_string());
        }
    }
}

impl std::fmt::Debug for RegionCache {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RegionCache")
            .field("len", &self.len())
            .field("capacity", &self.capacity)
            .finish()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_region(msg: &str) -> CachedRegion {
        CachedRegion {
            ref_data: CachedPanelData::default(),
            hap1_data: CachedPanelData::default(),
            hap2_data: CachedPanelData::default(),
            status_message: msg.to_string(),
        }
    }

    #[test]
    fn test_cache_new() {
        let cache = RegionCache::new(16);
        assert_eq!(cache.len(), 0);
        assert!(cache.is_empty());
        assert_eq!(cache.capacity(), 16);
    }

    #[test]
    fn test_cache_insert_and_get() {
        let mut cache = RegionCache::new(4);
        cache.insert("chr1:100-200".to_string(), dummy_region("first"));
        assert_eq!(cache.len(), 1);

        let entry = cache.get("chr1:100-200");
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().status_message, "first");
    }

    #[test]
    fn test_cache_miss() {
        let mut cache = RegionCache::new(4);
        cache.insert("chr1:100-200".to_string(), dummy_region("first"));
        assert!(cache.get("chr2:300-400").is_none());
    }

    #[test]
    fn test_cache_eviction_lru() {
        let mut cache = RegionCache::new(3);
        cache.insert("a".to_string(), dummy_region("A"));
        cache.insert("b".to_string(), dummy_region("B"));
        cache.insert("c".to_string(), dummy_region("C"));
        assert_eq!(cache.len(), 3);

        // Insert a 4th — should evict "a" (least recently used)
        cache.insert("d".to_string(), dummy_region("D"));
        assert_eq!(cache.len(), 3);
        assert!(cache.get("a").is_none(), "a should have been evicted");
        assert!(cache.get("b").is_some());
        assert!(cache.get("c").is_some());
        assert!(cache.get("d").is_some());
    }

    #[test]
    fn test_cache_access_promotes() {
        let mut cache = RegionCache::new(3);
        cache.insert("a".to_string(), dummy_region("A"));
        cache.insert("b".to_string(), dummy_region("B"));
        cache.insert("c".to_string(), dummy_region("C"));

        // Access "a" to promote it
        let _ = cache.get("a");

        // Insert "d" — should evict "b" now (it is the LRU after "a" was promoted)
        cache.insert("d".to_string(), dummy_region("D"));
        assert_eq!(cache.len(), 3);
        assert!(cache.get("a").is_some(), "a should be retained (promoted)");
        assert!(cache.get("b").is_none(), "b should have been evicted");
    }

    #[test]
    fn test_cache_update_existing() {
        let mut cache = RegionCache::new(4);
        cache.insert("a".to_string(), dummy_region("old"));
        cache.insert("a".to_string(), dummy_region("new"));
        assert_eq!(cache.len(), 1);
        assert_eq!(cache.get("a").unwrap().status_message, "new");
    }

    #[test]
    fn test_cache_clear() {
        let mut cache = RegionCache::new(4);
        cache.insert("a".to_string(), dummy_region("A"));
        cache.insert("b".to_string(), dummy_region("B"));
        assert_eq!(cache.len(), 2);
        cache.clear();
        assert_eq!(cache.len(), 0);
        assert!(cache.is_empty());
        assert!(cache.get("a").is_none());
    }

    #[test]
    fn test_cache_capacity_one() {
        let mut cache = RegionCache::new(1);
        cache.insert("a".to_string(), dummy_region("A"));
        cache.insert("b".to_string(), dummy_region("B"));
        assert_eq!(cache.len(), 1);
        assert!(cache.get("a").is_none());
        assert!(cache.get("b").is_some());
    }

    #[test]
    fn test_cache_debug_format() {
        let cache = RegionCache::new(8);
        let dbg = format!("{cache:?}");
        assert!(dbg.contains("RegionCache"));
        assert!(dbg.contains("capacity: 8"));
    }
}
