mod implementation;
pub use implementation::{Cluster, Hotspot, HotspotClass, find_hotspots, write_pdbstr};

mod python;
use pyo3::prelude::{Bound, PyModule, PyModuleMethods, PyResult, pymodule, wrap_pyfunction};

#[pymodule]
#[pyo3(name = "xdrugpy_hotspot_finder")]
pub fn register_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<python::PyHotspotClass>()?;
    m.add_class::<python::PyCluster>()?;
    m.add_class::<python::PyHotspot>()?;
    m.add_function(wrap_pyfunction!(python::py_find_hotspots, m)?)?;
    m.add_function(wrap_pyfunction!(python::py_write_pdbstr, m)?)?;
    Ok(())
}
