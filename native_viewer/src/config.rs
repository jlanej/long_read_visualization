use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};

use crate::genome::region::GenomicRegion;

/// A single sample's file paths and metadata.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleConfig {
    pub sample_id: String,
    pub output_dir: Option<PathBuf>,
    pub reference: Option<PathBuf>,
    pub hap1_assembly: Option<PathBuf>,
    pub hap2_assembly: Option<PathBuf>,
    pub reads_bam: Option<PathBuf>,
    pub regions_file: Option<PathBuf>,
    pub cram_ref: Option<PathBuf>,
    // Derived paths (auto-discovered from output_dir)
    pub hap1_to_ref_bam: Option<PathBuf>,
    pub hap2_to_ref_bam: Option<PathBuf>,
    pub reads_to_hap1_bam: Option<PathBuf>,
    pub reads_to_hap2_bam: Option<PathBuf>,
    pub hap1_mapping_index: Option<PathBuf>,
    pub hap2_mapping_index: Option<PathBuf>,
}

/// Manifest region entry (from toy_manifest.json or similar).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ManifestRegion {
    pub chrom: String,
    pub pos: u64,
    pub size: u64,
    #[serde(default)]
    pub genotype: Option<String>,
    #[serde(default)]
    pub ref_region: Option<String>,
    #[serde(default)]
    pub fasta_region: Option<String>,
    #[serde(default)]
    pub hap1_regions: Vec<String>,
    #[serde(default)]
    pub hap2_regions: Vec<String>,
}

/// Manifest JSON file format.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Manifest {
    #[serde(default)]
    pub description: Option<String>,
    pub variants: Vec<ManifestRegion>,
}

/// Application configuration encompassing all samples and regions.
#[derive(Debug, Clone)]
pub struct AppConfig {
    pub samples: Vec<SampleConfig>,
    pub regions: Vec<GenomicRegion>,
}

impl AppConfig {
    /// Load configuration from a TSV file (same format as the Python server).
    ///
    /// TSV columns (tab-separated, first line may start with `#`):
    /// sample_id, output_dir, reference, hap1_assembly, hap2_assembly,
    /// reads_bam, regions, cram_ref
    pub fn from_tsv(path: &str) -> Result<Self> {
        let content = fs::read_to_string(path)
            .with_context(|| format!("Failed to read config TSV: {path}"))?;

        let mut samples = Vec::new();
        for line in content.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 2 {
                continue;
            }

            let sample = SampleConfig {
                sample_id: cols[0].to_string(),
                output_dir: non_empty_path(cols.get(1)),
                reference: non_empty_path(cols.get(2)),
                hap1_assembly: non_empty_path(cols.get(3)),
                hap2_assembly: non_empty_path(cols.get(4)),
                reads_bam: non_empty_path(cols.get(5)),
                regions_file: non_empty_path(cols.get(6)),
                cram_ref: non_empty_path(cols.get(7)),
                hap1_to_ref_bam: None,
                hap2_to_ref_bam: None,
                reads_to_hap1_bam: None,
                reads_to_hap2_bam: None,
                hap1_mapping_index: None,
                hap2_mapping_index: None,
            };
            samples.push(sample);
        }

        // Auto-discover pipeline output files
        for sample in &mut samples {
            discover_pipeline_files(sample);
        }

        // Load regions from the first sample's regions file
        let regions = if let Some(sample) = samples.first() {
            load_regions(sample)?
        } else {
            Vec::new()
        };

        Ok(Self { samples, regions })
    }

    /// Create config from individual CLI arguments.
    pub fn from_cli(
        reference: Option<&str>,
        bam: Option<&str>,
        regions: Option<&str>,
    ) -> Self {
        let sample = SampleConfig {
            sample_id: "cli".to_string(),
            output_dir: None,
            reference: reference.map(PathBuf::from),
            hap1_assembly: None,
            hap2_assembly: None,
            reads_bam: bam.map(PathBuf::from),
            regions_file: regions.map(PathBuf::from),
            cram_ref: None,
            hap1_to_ref_bam: None,
            hap2_to_ref_bam: None,
            reads_to_hap1_bam: None,
            reads_to_hap2_bam: None,
            hap1_mapping_index: None,
            hap2_mapping_index: None,
        };

        let regions = load_regions(&sample).unwrap_or_default();

        Self {
            samples: vec![sample],
            regions,
        }
    }

    /// Number of configured regions of interest.
    pub fn region_count(&self) -> usize {
        self.regions.len()
    }

    /// Get a region by index.
    pub fn get_region(&self, idx: usize) -> Option<&GenomicRegion> {
        self.regions.get(idx)
    }
}

/// Auto-discover pipeline output files from a sample's output directory.
fn discover_pipeline_files(sample: &mut SampleConfig) {
    let Some(output_dir) = &sample.output_dir else {
        return;
    };

    if !output_dir.is_dir() {
        return;
    }

    // Look for hap1/hap2_to_ref BAMs and mapping indices
    if let Ok(entries) = fs::read_dir(output_dir) {
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().to_string();
            let path = entry.path();

            if name.ends_with("_hap1_to_ref.bam") {
                sample.hap1_to_ref_bam = Some(path);
            } else if name.ends_with("_hap2_to_ref.bam") {
                sample.hap2_to_ref_bam = Some(path);
            } else if name.ends_with("_reads_to_hap1.bam") || name.ends_with("_reads_to_hap1.cram")
            {
                sample.reads_to_hap1_bam = Some(path);
            } else if name.ends_with("_reads_to_hap2.bam") || name.ends_with("_reads_to_hap2.cram")
            {
                sample.reads_to_hap2_bam = Some(path);
            } else if name.ends_with("_hap1_to_ref.mapping.json.gz") {
                sample.hap1_mapping_index = Some(path);
            } else if name.ends_with("_hap2_to_ref.mapping.json.gz") {
                sample.hap2_mapping_index = Some(path);
            }
        }
    }
}

/// Load regions from a manifest JSON or VCF file.
fn load_regions(sample: &SampleConfig) -> Result<Vec<GenomicRegion>> {
    let Some(regions_file) = &sample.regions_file else {
        return Ok(Vec::new());
    };

    if !regions_file.exists() {
        return Ok(Vec::new());
    }

    let path_str = regions_file.to_string_lossy();

    if path_str.ends_with(".json") {
        load_manifest_regions(regions_file)
    } else if path_str.ends_with(".vcf") || path_str.ends_with(".vcf.gz") {
        // VCF support placeholder - would parse VCF for SV regions
        Ok(Vec::new())
    } else {
        Ok(Vec::new())
    }
}

/// Load regions from a manifest JSON file.
fn load_manifest_regions(path: &Path) -> Result<Vec<GenomicRegion>> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("Failed to read manifest: {}", path.display()))?;
    let manifest: Manifest =
        serde_json::from_str(&content).with_context(|| "Failed to parse manifest JSON")?;

    let regions: Vec<GenomicRegion> = manifest
        .variants
        .iter()
        .filter_map(|v| {
            // Use ref_region if available, otherwise construct from chrom/pos/size
            if let Some(ref region_str) = v.ref_region {
                GenomicRegion::parse(region_str)
            } else {
                Some(GenomicRegion {
                    chrom: v.chrom.clone(),
                    start: v.pos,
                    end: v.pos + v.size,
                })
            }
        })
        .collect();

    Ok(regions)
}

/// Convert an optional TSV column to a PathBuf if non-empty.
fn non_empty_path(col: Option<&&str>) -> Option<PathBuf> {
    col.and_then(|s| {
        let s = s.trim();
        if s.is_empty() { None } else { Some(PathBuf::from(s)) }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_manifest_region() {
        let json = r#"{
            "description": "test",
            "variants": [
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "size": 500,
                    "ref_region": "chr1:100-600",
                    "hap1_regions": ["hap1_contig:50-550"],
                    "hap2_regions": []
                }
            ]
        }"#;
        let manifest: Manifest = serde_json::from_str(json).unwrap();
        assert_eq!(manifest.variants.len(), 1);
        assert_eq!(manifest.variants[0].chrom, "chr1");
        assert_eq!(manifest.variants[0].pos, 100);
        assert_eq!(manifest.variants[0].size, 500);
    }

    #[test]
    fn test_non_empty_path() {
        assert!(non_empty_path(Some(&"")).is_none());
        assert!(non_empty_path(Some(&"  ")).is_none());
        assert!(non_empty_path(None).is_none());
        assert_eq!(
            non_empty_path(Some(&"/path/to/file")),
            Some(PathBuf::from("/path/to/file"))
        );
    }

    #[test]
    fn test_from_cli_empty() {
        let config = AppConfig::from_cli(None, None, None);
        assert_eq!(config.samples.len(), 1);
        assert_eq!(config.samples[0].sample_id, "cli");
        assert!(config.regions.is_empty());
    }

    #[test]
    fn test_load_manifest_regions() {
        let dir = tempfile::tempdir().unwrap();
        let manifest_path = dir.path().join("manifest.json");
        fs::write(
            &manifest_path,
            r#"{"variants":[
                {"chrom":"chr1","pos":1000,"size":500,"ref_region":"chr1:1000-1500"},
                {"chrom":"chr2","pos":2000,"size":300}
            ]}"#,
        )
        .unwrap();

        let regions = load_manifest_regions(&manifest_path).unwrap();
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].chrom, "chr1");
        assert_eq!(regions[0].start, 1000);
        assert_eq!(regions[0].end, 1500);
        assert_eq!(regions[1].chrom, "chr2");
        assert_eq!(regions[1].start, 2000);
        assert_eq!(regions[1].end, 2300);
    }
}
