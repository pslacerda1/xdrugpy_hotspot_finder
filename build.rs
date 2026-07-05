fn main() {
    let version = std::env::var("__VERSION__").unwrap_or("@@DEVELOPER_ONLY@@".to_string());
    println!("cargo:rustc-env=__VERSION__={version}");
}
