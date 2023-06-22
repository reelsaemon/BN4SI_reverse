# BN4SI
Implementation of Process-aware Bayesian Networks for Sequential Inference on Event Logs

## Framework
We provide a framework to generate Process-aware Bayesian Network structures from event log data. A given dataset with case identifiers, activities and an ordering of activities is transformed to fit the model assumptions of Bayesian Networks for the underlying model structures to be free of cycles.

For that, the activity information in the trace data is recoded to handle loop structures inside traces. Moreover, the ordering information of activities is extracted with use of Eventually Follows Graphs. The extracted ordering relations are then checked for transitivity to ensure that possible cycle structures stemming from concurrency across multiple different traces are decomposed into linear structures. For the generated structure, we preprocess the event log data accordingly to ensure that parameter learning for the designed approach is possible. Eventually, unknown trace information can be predicted via Bayesian Inference by querying the Process-aware Bayesian Network with a given set of evidence about known prefixes or suffixes of traces.

## Usage
- Load a desired event log data set. We provide an exemplary synthetic event log (`data/sim_event_log.csv`).
- Adjust the predefined model setup and constants if necessary.
- Generate the Process-aware Network structure with the provided set of functions.
- Preprocess the data with the provided set of functions.
- Query the Process-aware Bayesian Network with the framework provided by the `bnlearn` package we rely on in our implementations.

## Evaluation Scripts
We provide two evaluation scripts. One for Next Activity and Remaining Trace Prediction (`evaluation_suffix.R`) and one for Previous Activity and Preceding Trace Prediction (`evaluation_prefix.R`).
