use approx::assert_relative_eq;
use std::fs::read_to_string;
use std::path::PathBuf;

use xdrugpy_hotspot_finder::*;

#[test]
fn test_general_loading_8agl() {
    let pdb_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("8AGL.pdb");
    let pdb_str = read_to_string(pdb_path).expect("can't read file");
    let (_, clusters, hotspots) =
        find_hotspots(pdb_str, 0.15, 25, 0.5, false, 15, false).expect("must have");

    assert_eq!(8, clusters.len());
    assert_eq!(hotspots.len(), 1);

    let hs = hotspots.first().expect("must have");
    assert_eq!(69, hs.strength_total);
    assert_relative_eq!(20.49, hs.max_distance, epsilon = 0.01);
    //// this is really undefined for 8AGL because there is clashes between CS0 and CS6
    // assert_relative_eq!(07.48, hs.centroid_distance.unwrap(), epsilon=0.01);
}

#[test]
fn test_general_loading_8b7j() {
    let pdb_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("8B7J.pdb");
    let pdb_str = read_to_string(pdb_path).expect("can't read file");
    let (_, clusters, hotspots) =
        find_hotspots(pdb_str, 0.10, 25, 0.5, true, 15, true).expect("must have");

    assert_eq!(9, clusters.len());
    assert_eq!(2, hotspots.len());

    let hs = hotspots.get(1).expect("must have");
    assert_eq!(46, hs.strength_total);
    assert_eq!(14, hs.strength_0);
    assert_eq!(13, hs.strength_1.unwrap());
    assert_eq!(5, hs.strength_z.unwrap());
    assert_relative_eq!(10.54, hs.centroid_distance.unwrap(), epsilon = 0.01);
    assert_relative_eq!(15.86, hs.max_distance, epsilon = 0.01);
    assert_eq!(5, hs.clusters.len());
}

#[test]
fn test_general_loading_cf_2tpr() {
    let pdb_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("Cf_2TPR.pdb");
    let pdb_str = read_to_string(pdb_path).expect("can't read file");
    let (_, clusters, hotspots) =
        find_hotspots(pdb_str, 0.10, 25, 0.5, true, 15, true).expect("must have");

    assert_eq!(8, clusters.len());
    assert_eq!(8, hotspots.len());

    let hs = hotspots.first().expect("must have");
    assert_eq!(58, hs.strength_total);
    assert_eq!(14, hs.strength_0);
    assert_eq!(12, hs.strength_1.unwrap());
    assert_eq!(5, hs.strength_z.unwrap());
    assert_relative_eq!(12.026, hs.centroid_distance.unwrap(), epsilon = 0.001);
    assert_relative_eq!(27.312, hs.max_distance, epsilon = 0.001);
    assert_eq!(7, hs.clusters.len());
}
