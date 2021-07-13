# IPGQC
A working repository for Identical Protein Group content assessment tool.
## Principle
Given a full list of Identical Protein Group entries (IPGs) for a selected taxon and a query containing translated reading frames from a genome assembly, the tool a) hashes the reference file and populates the SQL database with the hashed entries, and b) hashes the query sequences 
## Use
1. Download *Main.cpp*, *md5.cpp*, and *md5.h* scripts. (The other scripts in the current repository attest to IPG content assessment and do not affect the tool's performance)
2. Compile *Main.cpp* using the default settings of the compiler installed on your machine.
3. Launch the compiled script as following: 

`./<scriptname> -d <database> -s <query>` - if script is run for the first time and no database is available 

`./<scriptname> -s <query>` - if the database is already built 

Here, `<database>` refers to the FASTA file for selected taxon(-a) and `<query>` refers to protein content of genome assembly in question. For sequential launches of IPGQC with the same database, only the first use needs *-d* parameter specified.
**NOTE**: The current version does not support differential database use. If you switch to another database, the previous one will be overriden. This is to be fixed soon.
## Auxiliary scripts
Script *IPG_stats_hists_upd.R* contains code used for figure production as presented at BiATA 2021 (Saint Petersburg, Russia, held online). The *scripts* directory contain Shell scripts used for genome assembly, respective IPG content and IPG statistics retrieval.
## Expected updates
The first major update is going to have the following features:
1. Main script - database construction independent of query search; multithreading; multiple database handling for independent launches;
2. Auxiliary scripts - a full suite for user-friendly IPG sequence and statistics retrieval.
## Reference
IPGQC has not been published. For project rationale, refer to the following abstract:
Yury V. Malovichko, Ruslan O. Alagov, Anton E. Shikov, Alexander V. Predeus, Anton A. Nizhnikov, Kirill S. Antonets. Identical Protein Group content and resequencing statistics as a naive metric of biological assembly quality: An evaluation study and draft tool implementation. From abstract book: Bioinformatics: From Algorithms to Applications (BiATA). July 12-13, Saint Petersburg, Russia.
