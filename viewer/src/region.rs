use std::fmt;
use std::io::Read;
use std::path::Path;

use serde::Deserialize;

/// Format a base-pair count as a human-readable string.
fn format_bp(bp: u64) -> String {
    if bp >= 1_000_000 {
        format!("{:.1} Mb", bp as f64 / 1_000_000.0)
    } else if bp >= 1000 {
        format!("{:.1} kb", bp as f64 / 1000.0)
    } else {
        format!("{bp} bp")
    }
}

/// A genomic region: chrom:start-end (1-based, inclusive).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl fmt::Display for GenomicRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
    }
}

impl GenomicRegion {
    /// Parse "chr1:100-200" into a GenomicRegion.
    pub fn parse(s: &str) -> Option<Self> {
        let (chrom, rest) = s.split_once(':')?;
        let (start_s, end_s) = rest.split_once('-')?;
        let start = start_s.replace(',', "").parse::<u64>().ok()?;
        let end = end_s.replace(',', "").parse::<u64>().ok()?;
        if start == 0 || end < start {
            return None;
        }
        Some(Self {
            chrom: chrom.to_string(),
            start,
            end,
        })
    }
}

/// A variant/region entry from a manifest JSON file.
#[derive(Debug, Clone, Deserialize)]
#[allow(dead_code)]
pub struct ManifestVariant {
    pub chrom: String,
    pub pos: u64,
    #[serde(default)]
    pub size: u64,
    #[serde(default)]
    pub genotype: String,
    #[serde(default)]
    pub ref_region: String,
    #[serde(default)]
    pub fasta_region: String,
    #[serde(default)]
    pub hap1_regions: Vec<String>,
    #[serde(default)]
    pub hap2_regions: Vec<String>,
}

/// Top-level manifest JSON structure.
#[derive(Debug, Clone, Deserialize)]
#[allow(dead_code)]
pub struct Manifest {
    #[serde(default)]
    pub description: String,
    pub variants: Vec<ManifestVariant>,
}

/// Summary of the region currently being displayed.
#[derive(Debug, Clone)]
pub struct RegionEntry {
    pub label: String,
    pub ref_region: GenomicRegion,
    pub hap1_region: Option<GenomicRegion>,
    pub hap2_region: Option<GenomicRegion>,
    /// Genotype string (e.g. "0|1", "1|1").
    pub genotype: String,
    /// SV size in base pairs.
    pub sv_size: u64,
    /// Free-form description from manifest (retained for future UI display).
    pub description: String,
}

impl fmt::Display for RegionEntry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.label)
    }
}

/// Manages a list of regions and a cursor for navigation.
#[derive(Debug)]
pub struct RegionNavigator {
    regions: Vec<RegionEntry>,
    current: usize,
}

impl Default for RegionNavigator {
    fn default() -> Self {
        Self::new()
    }
}

impl RegionNavigator {
    pub fn new() -> Self {
        Self {
            regions: Vec::new(),
            current: 0,
        }
    }

    /// Load regions from a file, dispatching by extension.
    ///
    /// Supported formats:
    /// - `.json` — manifest JSON (existing format)
    /// - `.vcf` — plain-text VCF
    /// - `.vcf.gz` — bgzip/gzip-compressed VCF
    pub fn load_regions(&mut self, path: &Path) -> Result<(), String> {
        let path_str = path.to_string_lossy();
        if path_str.ends_with(".vcf.gz") {
            let file = std::fs::File::open(path).map_err(|e| format!("Failed to open VCF: {e}"))?;
            let mut decoder = flate2::read::MultiGzDecoder::new(file);
            let mut text = String::new();
            decoder
                .read_to_string(&mut text)
                .map_err(|e| format!("Failed to decompress VCF: {e}"))?;
            self.load_vcf_str(&text)
        } else if path_str.ends_with(".vcf") {
            let text =
                std::fs::read_to_string(path).map_err(|e| format!("Failed to read VCF: {e}"))?;
            self.load_vcf_str(&text)
        } else {
            // Default: treat as manifest JSON
            let data = std::fs::read_to_string(path)
                .map_err(|e| format!("Failed to read manifest: {e}"))?;
            self.load_manifest_str(&data)
        }
    }

    /// Load regions from a manifest JSON file.
    #[cfg(test)]
    pub fn load_manifest(&mut self, path: &Path) -> Result<(), String> {
        let data =
            std::fs::read_to_string(path).map_err(|e| format!("Failed to read manifest: {e}"))?;
        self.load_manifest_str(&data)
    }

    /// Load regions from a manifest JSON string.
    pub fn load_manifest_str(&mut self, json: &str) -> Result<(), String> {
        let manifest: Manifest =
            serde_json::from_str(json).map_err(|e| format!("Invalid manifest JSON: {e}"))?;
        self.regions.clear();
        self.current = 0;

        for v in &manifest.variants {
            let ref_region = if !v.ref_region.is_empty() {
                GenomicRegion::parse(&v.ref_region)
            } else {
                None
            };

            let ref_region = match ref_region {
                Some(r) => r,
                None => {
                    // Fallback: construct from chrom+pos+size
                    if v.size > 0 {
                        GenomicRegion {
                            chrom: v.chrom.clone(),
                            start: v.pos,
                            end: v.pos + v.size,
                        }
                    } else {
                        continue; // Skip malformed entries
                    }
                }
            };

            let hap1_region = v.hap1_regions.first().and_then(|s| GenomicRegion::parse(s));
            let hap2_region = v.hap2_regions.first().and_then(|s| GenomicRegion::parse(s));

            let label = format!(
                "{}:{} ({}bp)",
                v.chrom,
                v.pos,
                if v.size > 0 {
                    v.size.to_string()
                } else {
                    "?".to_string()
                }
            );

            self.regions.push(RegionEntry {
                label,
                ref_region,
                hap1_region,
                hap2_region,
                genotype: v.genotype.clone(),
                sv_size: v.size,
                description: manifest.description.clone(),
            });
        }

        Ok(())
    }

    /// Load regions from VCF text content.
    ///
    /// Parses VCF lines for SV records. Determines SV size from the `SVLEN`
    /// INFO field, falling back to `abs(len(REF) - len(ALT))`. Records with
    /// size < 500 bp are skipped (matching the Python server behaviour).
    pub fn load_vcf_str(&mut self, vcf_text: &str) -> Result<(), String> {
        const MIN_SIZE: u64 = 500;

        self.regions.clear();
        self.current = 0;

        for line in vcf_text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 8 {
                continue;
            }

            let chrom = cols[0].to_string();
            let pos: u64 = match cols[1].parse() {
                Ok(v) => v,
                Err(_) => continue,
            };
            let ref_allele = cols[3];
            let alt_allele = cols[4];
            let info = cols[7];

            // Determine SV size: prefer SVLEN from INFO, then allele diff.
            let mut sv_len: u64 = 0;
            for field in info.split(';') {
                if let Some(val) = field.strip_prefix("SVLEN=")
                    && let Ok(v) = val.parse::<i64>()
                {
                    sv_len = v.unsigned_abs();
                    break;
                }
            }
            if sv_len == 0 {
                let ref_len = ref_allele.len() as u64;
                let alt_len = alt_allele.len() as u64;
                sv_len = ref_len.abs_diff(alt_len);
            }
            if sv_len < MIN_SIZE {
                continue;
            }

            let end = pos + sv_len;

            // Parse genotype if available (column 10).
            let gt = if cols.len() >= 10 {
                cols[9].split(':').next().unwrap_or("").replace('/', "|")
            } else {
                String::new()
            };

            let mut label = format!("{chrom}:{pos}-{end} ({}", format_bp(sv_len));
            if !gt.is_empty() {
                label.push_str(&format!(", {gt}"));
            }
            label.push(')');

            self.regions.push(RegionEntry {
                label,
                ref_region: GenomicRegion {
                    chrom,
                    start: pos,
                    end,
                },
                hap1_region: None,
                hap2_region: None,
                genotype: gt,
                sv_size: sv_len,
                description: String::new(),
            });
        }

        Ok(())
    }

    /// Number of loaded regions.
    pub fn len(&self) -> usize {
        self.regions.len()
    }

    /// Whether there are no regions loaded.
    pub fn is_empty(&self) -> bool {
        self.regions.is_empty()
    }

    /// Current 0-based index.
    pub fn current_index(&self) -> usize {
        self.current
    }

    /// Get the current region entry, if any.
    pub fn current(&self) -> Option<&RegionEntry> {
        self.regions.get(self.current)
    }

    /// Navigate to next region. Wraps around.
    pub fn next(&mut self) {
        if !self.regions.is_empty() {
            self.current = (self.current + 1) % self.regions.len();
        }
    }

    /// Navigate to previous region. Wraps around.
    pub fn prev(&mut self) {
        if !self.regions.is_empty() {
            if self.current == 0 {
                self.current = self.regions.len() - 1;
            } else {
                self.current -= 1;
            }
        }
    }

    /// Navigate to a specific index (clamped).
    pub fn go_to(&mut self, index: usize) {
        if !self.regions.is_empty() {
            self.current = index.min(self.regions.len() - 1);
        }
    }

    /// Get the list of region labels for a dropdown/selector.
    pub fn labels(&self) -> Vec<String> {
        self.regions.iter().map(|r| r.label.clone()).collect()
    }

    /// Human-readable counter string: "1/10".
    pub fn counter_text(&self) -> String {
        if self.regions.is_empty() {
            "0/0".to_string()
        } else {
            format!("{}/{}", self.current + 1, self.regions.len())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn sample_manifest_json() -> &'static str {
        r#"{
          "description": "Test manifest",
          "variants": [
            {
              "chrom": "chr1",
              "pos": 100000,
              "size": 5000,
              "genotype": "0|1",
              "ref_region": "chr1:100000-105000",
              "fasta_region": "chr1:90000-115000",
              "hap1_regions": ["sample#1#ctg1:50000-55000"],
              "hap2_regions": ["sample#2#ctg2:60000-65000"]
            },
            {
              "chrom": "chr2",
              "pos": 200000,
              "size": 3000,
              "genotype": "1|0",
              "ref_region": "chr2:200000-203000",
              "fasta_region": "chr2:190000-213000",
              "hap1_regions": ["sample#1#ctg3:70000-73000"],
              "hap2_regions": ["sample#2#ctg4:80000-83000"]
            },
            {
              "chrom": "chr3",
              "pos": 300000,
              "size": 1000,
              "genotype": "1|1",
              "ref_region": "chr3:300000-301000",
              "fasta_region": "chr3:295000-306000",
              "hap1_regions": [],
              "hap2_regions": []
            }
          ]
        }"#
    }

    #[test]
    fn test_genomic_region_parse_valid() {
        let r = GenomicRegion::parse("chr1:100-200").unwrap();
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, 100);
        assert_eq!(r.end, 200);
    }

    #[test]
    fn test_genomic_region_parse_with_commas() {
        let r = GenomicRegion::parse("chr1:1,000-2,000").unwrap();
        assert_eq!(r.start, 1000);
        assert_eq!(r.end, 2000);
    }

    #[test]
    fn test_genomic_region_parse_invalid() {
        assert!(GenomicRegion::parse("chr1").is_none());
        assert!(GenomicRegion::parse("chr1:abc-200").is_none());
        assert!(GenomicRegion::parse("chr1:200-100").is_none()); // end < start
        assert!(GenomicRegion::parse("chr1:0-100").is_none()); // 0 start
        assert!(GenomicRegion::parse("").is_none());
    }

    #[test]
    fn test_genomic_region_display() {
        let r = GenomicRegion {
            chrom: "chr1".into(),
            start: 100,
            end: 200,
        };
        assert_eq!(r.to_string(), "chr1:100-200");
    }

    #[test]
    fn test_navigator_empty() {
        let nav = RegionNavigator::new();
        assert!(nav.is_empty());
        assert_eq!(nav.len(), 0);
        assert_eq!(nav.current_index(), 0);
        assert!(nav.current().is_none());
        assert_eq!(nav.counter_text(), "0/0");
    }

    #[test]
    fn test_navigator_load_manifest_str() {
        let mut nav = RegionNavigator::new();
        nav.load_manifest_str(sample_manifest_json()).unwrap();
        assert_eq!(nav.len(), 3);
        assert_eq!(nav.current_index(), 0);
        assert_eq!(nav.counter_text(), "1/3");
    }

    #[test]
    fn test_navigator_next_wraps() {
        let mut nav = RegionNavigator::new();
        nav.load_manifest_str(sample_manifest_json()).unwrap();

        assert_eq!(nav.current_index(), 0);
        nav.next();
        assert_eq!(nav.current_index(), 1);
        nav.next();
        assert_eq!(nav.current_index(), 2);
        nav.next();
        assert_eq!(nav.current_index(), 0); // wrap
    }

    #[test]
    fn test_navigator_prev_wraps() {
        let mut nav = RegionNavigator::new();
        nav.load_manifest_str(sample_manifest_json()).unwrap();

        assert_eq!(nav.current_index(), 0);
        nav.prev();
        assert_eq!(nav.current_index(), 2); // wrap to end
        nav.prev();
        assert_eq!(nav.current_index(), 1);
    }

    #[test]
    fn test_navigator_go_to() {
        let mut nav = RegionNavigator::new();
        nav.load_manifest_str(sample_manifest_json()).unwrap();

        nav.go_to(2);
        assert_eq!(nav.current_index(), 2);
        nav.go_to(100); // clamped
        assert_eq!(nav.current_index(), 2);
        nav.go_to(0);
        assert_eq!(nav.current_index(), 0);
    }

    #[test]
    fn test_navigator_labels() {
        let mut nav = RegionNavigator::new();
        nav.load_manifest_str(sample_manifest_json()).unwrap();

        let labels = nav.labels();
        assert_eq!(labels.len(), 3);
        assert!(labels[0].contains("chr1:100000"));
        assert!(labels[1].contains("chr2:200000"));
        assert!(labels[2].contains("chr3:300000"));
    }

    #[test]
    fn test_navigator_current_region_data() {
        let mut nav = RegionNavigator::new();
        nav.load_manifest_str(sample_manifest_json()).unwrap();

        let entry = nav.current().unwrap();
        assert_eq!(entry.ref_region.chrom, "chr1");
        assert_eq!(entry.ref_region.start, 100000);
        assert_eq!(entry.ref_region.end, 105000);
        assert!(entry.hap1_region.is_some());
        assert!(entry.hap2_region.is_some());

        // Third entry has no hap regions
        nav.go_to(2);
        let entry = nav.current().unwrap();
        assert!(entry.hap1_region.is_none());
        assert!(entry.hap2_region.is_none());
    }

    #[test]
    fn test_navigator_next_prev_on_empty() {
        let mut nav = RegionNavigator::new();
        nav.next(); // should not panic
        nav.prev(); // should not panic
        assert_eq!(nav.current_index(), 0);
    }

    #[test]
    fn test_navigator_load_manifest_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("manifest.json");
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(sample_manifest_json().as_bytes()).unwrap();

        let mut nav = RegionNavigator::new();
        nav.load_manifest(&path).unwrap();
        assert_eq!(nav.len(), 3);
    }

    #[test]
    fn test_navigator_load_invalid_json() {
        let mut nav = RegionNavigator::new();
        let result = nav.load_manifest_str("not json");
        assert!(result.is_err());
    }

    #[test]
    fn test_navigator_load_replaces_previous() {
        let mut nav = RegionNavigator::new();
        nav.load_manifest_str(sample_manifest_json()).unwrap();
        assert_eq!(nav.len(), 3);
        nav.next();
        assert_eq!(nav.current_index(), 1);

        // Reload resets cursor
        nav.load_manifest_str(sample_manifest_json()).unwrap();
        assert_eq!(nav.current_index(), 0);
        assert_eq!(nav.len(), 3);
    }

    #[test]
    fn test_navigator_counter_text() {
        let mut nav = RegionNavigator::new();
        assert_eq!(nav.counter_text(), "0/0");

        nav.load_manifest_str(sample_manifest_json()).unwrap();
        assert_eq!(nav.counter_text(), "1/3");
        nav.next();
        assert_eq!(nav.counter_text(), "2/3");
        nav.next();
        assert_eq!(nav.counter_text(), "3/3");
    }

    #[test]
    fn test_manifest_missing_optional_fields() {
        let json = r#"{
          "variants": [
            {
              "chrom": "chr1",
              "pos": 500,
              "size": 100,
              "ref_region": "chr1:500-600"
            }
          ]
        }"#;
        let mut nav = RegionNavigator::new();
        nav.load_manifest_str(json).unwrap();
        assert_eq!(nav.len(), 1);
        let entry = nav.current().unwrap();
        assert!(entry.hap1_region.is_none());
        assert!(entry.hap2_region.is_none());
    }

    #[test]
    fn test_manifest_fallback_to_chrom_pos_size() {
        let json = r#"{
          "variants": [
            {
              "chrom": "chr5",
              "pos": 1000,
              "size": 200,
              "ref_region": ""
            }
          ]
        }"#;
        let mut nav = RegionNavigator::new();
        nav.load_manifest_str(json).unwrap();
        assert_eq!(nav.len(), 1);
        let entry = nav.current().unwrap();
        assert_eq!(entry.ref_region.chrom, "chr5");
        assert_eq!(entry.ref_region.start, 1000);
        assert_eq!(entry.ref_region.end, 1200);
    }

    #[test]
    fn test_hap_region_parsing_with_hash_contig() {
        // hap regions use format: sample#hapN#contigName:start-end
        let r = GenomicRegion::parse("sample#1#ctg1:50000-55000");
        assert!(r.is_some());
        let r = r.unwrap();
        assert_eq!(r.chrom, "sample#1#ctg1");
        assert_eq!(r.start, 50000);
        assert_eq!(r.end, 55000);
    }

    #[test]
    fn test_load_repo_toy_manifest() {
        // Load the actual toy_manifest.json shipped with the repository
        let manifest_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("resources/toy_dataset/toy_manifest.json");
        if !manifest_path.exists() {
            // Skip if running outside repo tree
            return;
        }
        let mut nav = RegionNavigator::new();
        nav.load_manifest(&manifest_path).unwrap();
        assert_eq!(nav.len(), 10, "toy_manifest.json should have 10 variants");

        // Verify we can navigate through all 10
        for i in 0..10 {
            assert_eq!(nav.current_index(), i);
            let entry = nav.current().unwrap();
            assert!(!entry.ref_region.chrom.is_empty());
            assert!(entry.ref_region.start > 0);
            assert!(entry.ref_region.end >= entry.ref_region.start);
            nav.next();
        }
        // After 10 next() calls, we should wrap back to 0
        assert_eq!(nav.current_index(), 0);
    }

    // -- VCF loading tests --

    fn sample_vcf() -> String {
        [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
            &format!(
                "chr1\t460749\t.\t{}\tA\t.\t.\tSVLEN=834\tGT\t0|1",
                "A".repeat(835)
            ),
            &format!("chr2\t100000\t.\t{}\tC\t.\t.\t.\tGT\t1/0", "C".repeat(1501)),
            // Small SV (< 500 bp) – should be skipped
            "chr3\t200000\t.\tATCG\tA\t.\t.\tSVLEN=3\tGT\t0|1",
        ]
        .join("\n")
    }

    #[test]
    fn test_load_vcf_str_basic() {
        let mut nav = RegionNavigator::new();
        nav.load_vcf_str(&sample_vcf()).unwrap();
        // Two qualifying SVs (the 3-bp one is skipped)
        assert_eq!(nav.len(), 2);

        let e0 = nav.current().unwrap();
        assert_eq!(e0.ref_region.chrom, "chr1");
        assert_eq!(e0.ref_region.start, 460749);
        assert_eq!(e0.ref_region.end, 460749 + 834);
        assert!(e0.label.contains("0|1"));
        assert!(e0.hap1_region.is_none());
        assert!(e0.hap2_region.is_none());

        nav.next();
        let e1 = nav.current().unwrap();
        assert_eq!(e1.ref_region.chrom, "chr2");
        assert_eq!(e1.ref_region.start, 100000);
        assert_eq!(e1.ref_region.end, 100000 + 1500);
        // Genotype "/" should be normalised to "|"
        assert!(e1.label.contains("1|0"));
    }

    #[test]
    fn test_load_vcf_str_svlen_priority() {
        // SVLEN takes priority over allele length diff
        let vcf = "chr1\t1000\t.\tA\tA\t.\t.\tSVLEN=5000\tGT\t1|1\n";
        let mut nav = RegionNavigator::new();
        nav.load_vcf_str(vcf).unwrap();
        assert_eq!(nav.len(), 1);
        let e = nav.current().unwrap();
        assert_eq!(e.ref_region.end, 1000 + 5000);
    }

    #[test]
    fn test_load_vcf_str_empty() {
        let mut nav = RegionNavigator::new();
        nav.load_vcf_str("").unwrap();
        assert!(nav.is_empty());
    }

    #[test]
    fn test_load_vcf_str_only_comments() {
        let mut nav = RegionNavigator::new();
        nav.load_vcf_str("##fileformat=VCFv4.2\n#CHROM\tPOS\n")
            .unwrap();
        assert!(nav.is_empty());
    }

    #[test]
    fn test_load_regions_json() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("manifest.json");
        std::fs::write(&path, sample_manifest_json()).unwrap();

        let mut nav = RegionNavigator::new();
        nav.load_regions(&path).unwrap();
        assert_eq!(nav.len(), 3);
    }

    #[test]
    fn test_load_regions_vcf() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.vcf");
        std::fs::write(&path, sample_vcf()).unwrap();

        let mut nav = RegionNavigator::new();
        nav.load_regions(&path).unwrap();
        assert_eq!(nav.len(), 2);
    }

    #[test]
    fn test_load_regions_vcf_gz() {
        use flate2::Compression;
        use flate2::write::GzEncoder;

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.vcf.gz");
        let file = std::fs::File::create(&path).unwrap();
        let mut enc = GzEncoder::new(file, Compression::default());
        enc.write_all(sample_vcf().as_bytes()).unwrap();
        enc.finish().unwrap();

        let mut nav = RegionNavigator::new();
        nav.load_regions(&path).unwrap();
        assert_eq!(nav.len(), 2);
    }

    #[test]
    fn test_format_bp_sizes() {
        assert_eq!(format_bp(100), "100 bp");
        assert_eq!(format_bp(999), "999 bp");
        assert_eq!(format_bp(1000), "1.0 kb");
        assert_eq!(format_bp(1500), "1.5 kb");
        assert_eq!(format_bp(999_999), "1000.0 kb");
        assert_eq!(format_bp(1_000_000), "1.0 Mb");
        assert_eq!(format_bp(2_500_000), "2.5 Mb");
    }

    // -- Variant metadata tests --

    #[test]
    fn test_manifest_metadata_genotype() {
        let mut nav = RegionNavigator::new();
        nav.load_manifest_str(sample_manifest_json()).unwrap();

        let e0 = nav.current().unwrap();
        assert_eq!(e0.genotype, "0|1");
        assert_eq!(e0.sv_size, 5000);
        assert_eq!(e0.description, "Test manifest");

        nav.next();
        let e1 = nav.current().unwrap();
        assert_eq!(e1.genotype, "1|0");
        assert_eq!(e1.sv_size, 3000);
    }

    #[test]
    fn test_vcf_metadata_genotype() {
        let mut nav = RegionNavigator::new();
        nav.load_vcf_str(&sample_vcf()).unwrap();

        let e0 = nav.current().unwrap();
        assert_eq!(e0.genotype, "0|1");
        assert_eq!(e0.sv_size, 834);

        nav.next();
        let e1 = nav.current().unwrap();
        assert_eq!(e1.genotype, "1|0");
        assert_eq!(e1.sv_size, 1500);
    }

    #[test]
    fn test_region_entry_default_metadata() {
        let entry = RegionEntry {
            label: "test".to_string(),
            ref_region: GenomicRegion {
                chrom: "chr1".into(),
                start: 1,
                end: 100,
            },
            hap1_region: None,
            hap2_region: None,
            genotype: String::new(),
            sv_size: 0,
            description: String::new(),
        };
        assert!(entry.genotype.is_empty());
        assert_eq!(entry.sv_size, 0);
        assert!(entry.description.is_empty());
    }
}
