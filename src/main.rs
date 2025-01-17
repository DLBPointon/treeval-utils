use human_panic::setup_panic;
use treeval_utils::run;

// https://doc.rust-lang.org/book/ch12-03-improving-error-handling-and-modularity.html#separation-of-concerns-for-binary-projects
fn main() {
    //  https://rust-cli.github.io/book/in-depth/human-communcation.html
    setup_panic!();
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    } else {
        println!("Done!");
    }
}
