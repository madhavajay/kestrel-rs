use kestrel::cli::{CliCommand, help_topic, parse_command};
use kestrel::runner::RunnerError;

fn main() -> anyhow::Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .try_init()
        .ok();

    match parse_command(std::env::args())? {
        CliCommand::Help(topic) => {
            print!("{}", help_topic(&topic));
            Ok(())
        }
        CliCommand::Run(runner) => match runner.run() {
            Err(RunnerError::PipelineNotImplemented) => {
                anyhow::bail!("the Rust Kestrel runner pipeline is not implemented yet")
            }
            result => Ok(result?),
        },
    }
}
