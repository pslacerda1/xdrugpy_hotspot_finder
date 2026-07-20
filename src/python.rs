//
// Gerado por I.A..
//     (& ajustado manualmente)
//

use crate::{Cluster, Hotspot, HotspotClass, find_hotspots, write_pdbstr};
use pyo3::prelude::*;

// ─── HotspotClass ────────────────────────────────────────────────────────────

#[pyclass(name = "HotspotClass", eq, eq_int, skip_from_py_object)]
#[derive(Clone, PartialEq)]
pub enum PyHotspotClass {
    D,
    DS,
    DL,
    B,
    BS,
    BL,
}

impl From<HotspotClass> for PyHotspotClass {
    fn from(c: HotspotClass) -> Self {
        match c {
            HotspotClass::D => PyHotspotClass::D,
            HotspotClass::DS => PyHotspotClass::DS,
            HotspotClass::DL => PyHotspotClass::DL,
            HotspotClass::B => PyHotspotClass::B,
            HotspotClass::BS => PyHotspotClass::BS,
            HotspotClass::BL => PyHotspotClass::BL,
        }
    }
}

// ─── Cluster ─────────────────────────────────────────────────────────────────

#[pyclass(name = "Cluster", skip_from_py_object)]
#[derive(Clone)]
pub struct PyCluster(pub(crate) Cluster);

#[pymethods]
impl PyCluster {
    #[getter]
    pub fn title(&self) -> String {
        self.0.title.clone()
    }
    #[getter]
    pub fn strength(&self) -> u32 {
        self.0.strength
    }

    pub fn get_pdbstr(&mut self, atom_offset: usize) -> String {
        self.0.get_pdbstr(atom_offset)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "Cluster(title={:?}, strength={})",
            self.0.title, self.0.strength
        )
    }
}

// ─── Hotspot ─────────────────────────────────────────────────────────────────

#[pyclass(name = "Hotspot", skip_from_py_object)]
#[derive(Clone)]
pub struct PyHotspot(pub(crate) Hotspot);

#[pymethods]
impl PyHotspot {
    #[getter]
    pub fn class(&self) -> PyHotspotClass {
        self.0.class.clone().into()
    }
    #[getter]
    pub fn strength_total(&self) -> u32 {
        self.0.strength_total
    }
    #[getter]
    pub fn strength_0(&self) -> u32 {
        self.0.strength_0
    }
    #[getter]
    pub fn strength_1(&self) -> Option<u32> {
        self.0.strength_1
    }
    #[getter]
    pub fn strength_z(&self) -> Option<u32> {
        self.0.strength_z
    }
    #[getter]
    pub fn max_distance(&self) -> f32 {
        self.0.max_distance
    }
    #[getter]
    pub fn centroid_distance(&self) -> Option<f32> {
        self.0.centroid_distance
    }

    #[getter]
    pub fn clusters(&self) -> Vec<PyCluster> {
        self.0
            .clusters
            .iter()
            .map(|c| PyCluster(c.clone()))
            .collect()
    }

    pub fn get_pdbstr(&mut self) -> String {
        self.0.get_pdbstr()
    }

    pub fn __repr__(&self) -> String {
        format!(
            "Hotspot(class={:?}, strength_total={}, clusters={})",
            self.0.class,
            self.0.strength_total,
            self.0.clusters.len()
        )
    }
}

// ─── free functions ───────────────────────────────────────────────────────────

#[pyfunction]
#[pyo3(name = "find_hotspots")]
pub fn py_find_hotspots(
    pdb_str: String,
    clash_threshold: f32,
    num_pseudoatoms: u32,
    pseudoatom_radius: f32,
    deep_search: bool,
    max_size: u32,
    remove_nested: bool,
) -> PyResult<(Vec<String>, Vec<PyCluster>, Vec<PyHotspot>)> {
    let (protein_lines, clusters, hotspots) = find_hotspots(
        pdb_str,
        clash_threshold,
        num_pseudoatoms,
        pseudoatom_radius,
        deep_search,
        max_size,
        remove_nested,
    )
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    let py_clusters = clusters.iter().map(|c| PyCluster(c.clone())).collect();
    let py_hotspots = hotspots.iter().map(|h| PyHotspot(h.clone())).collect();
    Ok((protein_lines, py_clusters, py_hotspots))
}

#[pyfunction]
#[pyo3(name = "write_pdbstr")]
pub fn py_write_pdbstr(
    group: &str,
    path: String,
    protein_lines: Vec<String>,
    clusters: Vec<PyRef<PyCluster>>,
    hotspots: Vec<PyRef<PyHotspot>>,
) -> PyResult<()> {
    let mut file = std::fs::File::create(&path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let raw_clusters: Vec<Cluster> = clusters.iter().map(|c| c.0.clone()).collect();
    let raw_hotspots: Vec<Hotspot> = hotspots.iter().map(|h| h.0.clone()).collect();

    write_pdbstr(group, &mut file, protein_lines, raw_clusters, raw_hotspots)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(())
}
