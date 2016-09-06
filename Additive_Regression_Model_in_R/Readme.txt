The Computer Activity databases are a collection of computer systems activity measures. The datawas collected from a Sun Sparcstation 20/712 with 128 Mbytes of memory running in a multi-useruniversity department. Users would typically be doing a large variety of tasks ranging fromaccessing the internet, editing files or running heavy cpu-oriented programs. The data wascollected continuously which was system activity gathered every 5 seconds.

System measures that form feature vectors are as follows:

1. lread - Reads (transfers per second ) between system memory and user memory.2. lwrite - writes (transfers per second) between system memory and user memory.3. scall - Number of system calls of all types per second.4. sread - Number of system read calls per second.5. swrite - Number of system write calls per second .6. fork - Number of system fork calls per second.7. exec - Number of system exec calls per second.8. rchar - Number of characters transferred per second by system read calls.9. wchar - Number of characters transfreed per second by system write calls.10. runqsz - Process run queue size.11. freemem - Number of memory pages available to user processes.12. freeswap - Number of disk blocks available for page swapping.13. usr - Portion of time (%) that cpus run in user mode.

The goal of regression/ function estimation is to predict the attribute usr i.e. the portion oftime that CPU runs in user mode using attributes 1-12.

R-Code_and_Results.pdf contains implementation of Generalized Linear regression model in R programming language for solving this regression problem. Preprocessing of the data is also carried out.

