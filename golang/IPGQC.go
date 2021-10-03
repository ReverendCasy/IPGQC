package main

import (
	"bufio"
	"crypto/md5"
	"database/sql"
	"encoding/hex"
	"flag"
	"fmt"
	_ "github.com/mattn/go-sqlite3"
	"github.com/shenwei356/bio/seqio/fastx"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

type Record struct {
	id   int
	hash string
}

type SeqHash struct {
	id   string
	hash string
}

type  SearchResult struct {
	ipgid int
	specie string
}

type ReturningResult struct {
	proteinsPassed int
	proteinsDB int
	species string
}

type SpeciesCount struct {
	name string
	qty int
}

type Records []Record

//TODO: add cmd parameters
var db = flag.String("dbfile", "data.db", "path to database file")
var proteins = flag.String("proteins", "test_seqs.fasta", "path to file with proteins")
var ipgidsProtein = flag.String("ipgid", "test_stats.tsv", "path to file with matching ipgid and sequence name in fasta")
var speciesIpgid = flag.String("species", "patable.tsv", "path to file with species and information about proteins they contain")
//var proteinsCheck = flag.String("search", "check_seqs.fasta", "path to file with searching sequences")
var proteinsCheck = flag.String("search", "sen_flye_SRP250949_pilon.faa", "path to file with searching sequences")

func main() {
	//fmt.Println(time.Now().Format("2006-01-02 15:04:05"))
	flag.Parse()
	dbName := *db
	proteinsFile := *proteins
	ipgidstatsFile := *ipgidsProtein
	speciesIpgidFile := *speciesIpgid
	proteinsCheckFile := *proteinsCheck
	db := InitDB(dbName)
	defer db.Close()
	createTables(db)
	if checkProteinNotExists(db) {
		//stats_file := "test_stats.tsv"
		stats_data := parseTSV(ipgidstatsFile)
		//file := "test_seqs.fasta"
		//outfh, err := xopen.Wopen(file)
		//checkError(err)
		//defer outfh.Close()
		hashes, err := readFastaFile(proteinsFile)
		checkError(err)
		records := prepareRecords(stats_data, hashes)
		writeToDB(db, records)
	}
	if checkSpeciesNotExists(db) {
		//species := "patable.tsv"
		species_data := parseTSV(speciesIpgidFile)
		writeSpecies(db, species_data)
	}
	searchingProteinsHashes, err := readFastaFile(proteinsCheckFile)
	//_ = proteinsCheckFile
	//searchingProteinsHashes, err := readFastaFile(proteinsFile)
	checkError(err)
	var res ReturningResult
	species, proteinsPassed := searchProtein(db, searchingProteinsHashes)
	res.proteinsPassed = proteinsPassed
	//species := searchSpeciesWithProtein(db, protein)
	//res.proteinsPassed = cap(searchingProteinsHashes) / 2
	var proteinsDB []int
	var speciesCount []SpeciesCount
	for _, specie := range species{
		if !containsInt(proteinsDB, specie.ipgid) {
			proteinsDB = append(proteinsDB, specie.ipgid)
		}
		if !containsSpecie(speciesCount, specie.specie) {
			var specieNew SpeciesCount
			specieNew.qty = 1
			specieNew.name = specie.specie
			speciesCount = append(speciesCount, specieNew)
		} else {
			for _, s := range speciesCount {
				if s.name == specie.specie {
					s.qty = s.qty + 1
				}
			}
		}
	}
	sort.SliceStable(speciesCount, func(i, j int) bool {
		return speciesCount[i].qty > speciesCount[j].qty
	})
	res.proteinsDB = len(proteinsDB)
	fmt.Println("Total proteins " + strconv.Itoa(res.proteinsPassed) + ". Proteins in db " + strconv.Itoa(res.proteinsDB) + ".")
	fmt.Println("Species:")
	for _, s := range speciesCount {
		fmt.Println(s.name)
	}
	//fmt.Println(time.Now().Format("2006-01-02 15:04:05"))
	//fmt.Println("Finish")
}

func containsString(s []string, searchterm string) bool {
	i := sort.SearchStrings(s, searchterm)
	return i < len(s) && s[i] == searchterm
}

func containsSpecie(species []SpeciesCount, specie string) bool {
	for _, s := range species{
		if s.name == specie {
			return true
		}
	}
	return false
}

func containsInt(s []int, e int) bool {
	for _, a := range s {
		if a == e {
			return true
		}
	}
	return false
}

func searchProtein(db *sql.DB, hashes <-chan *SeqHash) ([]SearchResult, int) {
	var res []SearchResult
	var searchingProteinString string
	var seqHashCount int = 0
	for hash := range hashes {
		searchingProteinString = searchingProteinString + "'" + hash.hash + "'" + ", "
		seqHashCount = seqHashCount + 1
	}
	searchingProteinString = strings.TrimSpace(searchingProteinString)
	searchingProteinString = strings.TrimSuffix(searchingProteinString, ",")
	sqlStatement := "select phd.id, ss.name " +
		"from ProteinHashed as phd " +
		"inner join ProteinHashedSpecies as phs on phd.id = phs.ProteinHashed " +
		"inner join Species as ss on phs.Species = ss.id " +
		"where phd.Hash in (" + searchingProteinString + ")"

	stmt, err := db.Prepare(sqlStatement)
	rows, err := stmt.Query()
	defer rows.Close()
	if err != nil {
		log.Fatalln(err)
	}
	for rows.Next() {
		var ipgid int
		var name string
		err = rows.Scan(&ipgid, &name)
		if err != nil {
			log.Fatalln(err)
		}
		var r SearchResult
		r.ipgid = ipgid
		r.specie = name
		res = append(res, r)
	}
	return res, seqHashCount
}

func searchSpeciesWithProtein(db *sql.DB, protein string) []string{
	sqlStatement := "select ss.name from ProteinHashed as phd inner join ProteinHashedSpecies as phs on phd.id = phs.ProteinHashed inner join Species as ss on phs.Species = ss.id where phd.Hash = $1"
	rows, err := db.Query(sqlStatement, protein)
	checkError(err)
	defer rows.Close()
	var species []string
	for rows.Next() {
		var name string
		err = rows.Scan(&name)
		if err != nil {
			log.Fatalln(err)
		}
		species = append(species, name)
	}
	return species
}

func checkProteinNotExists(db *sql.DB) bool {
	sqlStatement := "select count(id) as proteinsCount from ProteinHashed"
	stmt, err := db.Prepare(sqlStatement)
	checkError(err)
	defer stmt.Close()
	tx, err := db.Begin()
	checkError(err)
	var proteinsCount int
	err = tx.Stmt(stmt).QueryRow().Scan(&proteinsCount)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		tx.Rollback()
		os.Exit(1)
	}
	tx.Commit()
	if proteinsCount == 0 {
		return true
	}
	return false
}

func checkSpeciesNotExists(db *sql.DB) bool {
	sqlStatement := "select count(id) as speciesCount from Species"
	stmt, err := db.Prepare(sqlStatement)
	checkError(err)
	defer stmt.Close()
	tx, err := db.Begin()
	checkError(err)
	var speciesCount int
	err = tx.Stmt(stmt).QueryRow().Scan(&speciesCount)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		tx.Rollback()
		os.Exit(1)
	}
	tx.Commit()
	if speciesCount == 0 {
		return true
	}
	return false
}

func readFastaFile(filename string) (<-chan *SeqHash, error) {
	res := make(chan *SeqHash, 10)
	reader, err := fastx.NewDefaultReader(filename)
	if err != nil {
		return res, err
	}
	go func() {
		defer reader.Close()
		for {
			record, err := reader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				log.Fatalln(err)
			}

			id := string(record.ID)
			id = id[:len(id)-2]
			hash := getMD5Hash(string(record.Seq.Seq))
			checkError(err)
			res <- &SeqHash{id, hash}
		}
		close(res)
	}()
	return res, err
}

func prepareRecords(stats_data map[string][]string, hashes <-chan *SeqHash) <-chan *Record {
	res := make(chan *Record, 10)
	go func() {
		for hash := range hashes {
			ipgId := stats_data[hash.id]
			s, err := strconv.Atoi(ipgId[1])
			if err != nil {
				log.Fatalln(err)
			}
			res <- &Record{
				id:   s,
				hash: hash.hash,
			}
		}
		close(res)
	}()
	return res
}

func checkError(err error) {
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func InitDB(dbName string) *sql.DB {
	if _, err := os.Stat(dbName); os.IsNotExist(err) {
		os.Create(dbName)
	}
	db, err := sql.Open("sqlite3", dbName)
	checkError(err)
	return db
}

func createTables(db *sql.DB) {
	var tableExists int
	row := db.QueryRow("select count(*) as tableExists from sqlite_master where type='table' AND name='ProteinHashed'")
	row.Scan(&tableExists)
	if tableExists == 0 {
		_, err := db.Exec("CREATE TABLE 'ProteinHashed' ('id' INTEGER PRIMARY KEY, 'Hash' CHAR (100))")
		checkError(err)
	}
	row = db.QueryRow("select count(*) as tableExists from sqlite_master where type='table' and name='Species'")
	row.Scan(&tableExists)
	if tableExists == 0 {
		_, err := db.Exec("create table 'Species' ('id' integer primary key, 'Name' char(200))")
		checkError(err)
	}
	row = db.QueryRow("select count(*) as tableExists from sqlite_master where type='table' and name='ProteinHashedSpecies'")
	row.Scan(&tableExists)
	if tableExists == 0 {
		_, err := db.Exec("create table 'ProteinHashedSpecies' ('ProteinHashed' integer, 'Species' integer, foreign key(ProteinHashed) references ProteinHashed(id), foreign key(Species) references Species(id))")
		checkError(err)
	}
}

func writeToDB(db *sql.DB, items <-chan *Record) {
	sqlStatement := "insert into ProteinHashed (id, Hash) values ($1, $2)"
	stmt, err := db.Prepare(sqlStatement)
	checkError(err)
	defer stmt.Close()
	tx, err := db.Begin()
	checkError(err)
	for item := range items {
		_, err := tx.Stmt(stmt).Exec(item.id, item.hash)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			tx.Rollback()
			os.Exit(1)
		}
	}
	tx.Commit()
}

func writeSpecies(db *sql.DB, species_records map[string][]string)  {
	sqlStatement := "insert into Species (Name) values ($1)"
	stmt, err := db.Prepare(sqlStatement)
	checkError(err)
	defer stmt.Close()
	sqlStatementIpgSpecie := "insert into ProteinHashedSpecies (ProteinHashed, Species) values ($1, $2)"
	stmtIpgS, err := db.Prepare(sqlStatementIpgSpecie)
	checkError(err)
	defer stmtIpgS.Close()
	tx, err := db.Begin()
	checkError(err)
	for item, specieIpgsids := range species_records {
		res, err := tx.Stmt(stmt).Exec(item)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			tx.Rollback()
			os.Exit(1)
		}
		specieId, err := res.LastInsertId()
		ipgids := strings.Split(specieIpgsids[1], ",")
		for _, ipgid := range ipgids {
			if ipgid != "" {
				id, err := strconv.Atoi(ipgid)
				if err != nil {
					fmt.Fprintln(os.Stderr, err)
					tx.Rollback()
					os.Exit(1)
				}
				_, err = tx.Stmt(stmtIpgS).Exec(id, specieId)
				if err != nil {
					fmt.Fprintln(os.Stderr, err)
					tx.Rollback()
					os.Exit(1)
				}
			}
		}
	}
	tx.Commit()
}

func parseTSV(fileName string) map[string][]string {
	tsv, err := os.Open(fileName)
	checkError(err)

	scn := bufio.NewScanner(tsv)
	const maxCapacity = 20000000
	buf := make([]byte, maxCapacity)
	scn.Buffer(buf, maxCapacity)

	var lines []string

	for scn.Scan() {
		line := scn.Text()
		lines = append(lines, line)
	}

	err = scn.Err()
	checkError(err)

	//lines = lines[1:] // First line is header
	out := make(map[string][]string, len(lines))

	for _, line := range lines {
		record := strings.Split(line, "\t")
		out[record[0]] = record
	}

	return out
}

func getMD5Hash(text string) string {
	hash := md5.Sum([]byte(text))
	return hex.EncodeToString(hash[:])
}
