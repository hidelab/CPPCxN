# PCxN_sharc_1/2
## Differences from the original Improved PCxN

1. Contains the pieces job array
2. Contains the parts job array
3. Slight changes in the PCxN_estimatesXX.R scripts mainly to accept the new $SGE_TASK_ID argument

NOTE: The parts job array was implemented by setting a very large number of jobs as the array's upper limit(e.g 300 or 500). This ensures that we will get execute every part even our genesets were too many and the script generated hundreds of parts to calculate. In the unlikely event that we have a gigantic number of parts (e.g. 600-700) just re-run the script with a job array range that will cover it.

Example:
Our genesets' file is so big that the script decided it needs 630 parts. We already ran the default PCxN_sharc_2 script which has range: `#$ -t 1-500`. It finished but we don't have 630 part.RDS files (naturally). We adjust the job array flag to `#$ -t 501-630` and run PCxN_sharc_2 again.
