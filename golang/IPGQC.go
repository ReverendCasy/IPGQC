package main

import (
	"bufio"
	"crypto/md5"
	"database/sql"
	"encoding/hex"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"

	fisher "github.com/glycerine/golang-fisher-exact"
	_ "github.com/mattn/go-sqlite3"
	"github.com/shenwei356/bio/seqio/fastx"
)

type Record struct {
	id   int
	hash string
}

type SeqHash struct {
	id   string
	hash string
}

type SearchResult struct {
	ipgid  int
	specie string
}

type ReturningResult struct {
	proteinsPassed int
	proteinsDB     int
	species        string
}

type SpeciesCount struct {
	name string
	qty  int
}

type Records []Record

var db = flag.String("dbfile", "data.db", "path to database file")
var proteinsFilePath = flag.String("proteinsFilePath", "test_seqs.fasta", "path to file with proteinsFilePath")
var ipgidsProtein = flag.String("ipgid", "test_stats.tsv", "path to file with matching ipgid and sequence name in fasta")
var speciesIpgid = flag.String("species", "patable.tsv", "path to file with species and information about proteinsFilePath they contain")

// var proteinsCheck = flag.String("search", "check_seqs.fasta", "path to file with searching sequences")
var proteinsCheck = flag.String("search", "sen_flye_SRP250949_pilon.faa", "path to file with searching sequences")
var onlyDB = flag.Bool("onlyDB", true, "create only database without search")

func main() {
	flag.Parse()
	dbName := *db
	proteinsFile := *proteinsFilePath
	ipgidstatsFile := *ipgidsProtein
	speciesIpgidFile := *speciesIpgid
	proteinsCheckFile := *proteinsCheck
	createOnlyDB := *onlyDB
	db := InitDB(dbName)
	defer db.Close() // the block with creating table should be refactored into another function and run as a special mode
	if createOnlyDB {
		createTables(db)
		if checkProteinNotExists(db) {
			stats_data := parseTSV(ipgidstatsFile)
			hashes, err := readFastaFile(proteinsFile)
			checkError(err)
			records := prepareRecords(stats_data, hashes)
			writeToDB(db, records)
		}
		if checkSpeciesNotExists(db) {
			species_data := parseTSV(speciesIpgidFile)
			writeSpecies(db, species_data)
		}
	}

	searchingProteinsHashes, err := readFastaFile(proteinsCheckFile)
	checkError(err)
	var totalProteins = len(searchingProteinsHashes)

	var res ReturningResult // the same. The block with per species stats should be moved to another function.
	proteinsDB, speciesQtyProtein := searchProtein(db, searchingProteinsHashes)
	res.proteinsPassed = len(searchingProteinsHashes)

	sort.SliceStable(speciesQtyProtein, func(i, j int) bool {
		return speciesQtyProtein[i].qty > speciesQtyProtein[j].qty
	})

	res.proteinsDB = len(proteinsDB)
	fmt.Println("Total proteins " + strconv.Itoa(res.proteinsPassed) + ".")
	fmt.Println("Proteins in db " + strconv.Itoa(res.proteinsDB) + ":")
	proteinsStr := ""
	for _, s := range proteinsDB {
		proteinsStr = proteinsStr + strconv.Itoa(s) + ", "
	}

	proteinsStr = strings.TrimSpace(proteinsStr)
	proteinsStr = strings.TrimSuffix(proteinsStr, ",")
	fmt.Println(proteinsStr)

	fmt.Println("Species:")
	for _, s := range speciesQtyProtein {
		fmt.Println(s.name + "	(" + strconv.Itoa(s.qty) + ")")
	}
	
	var totalDBProteins = countDBProteins(db)
	m := totalProteins - res.proteinsPassed                              // m should be proteinPassed - proteins attributed to the speciesQtyProtein
	var expSpec = (totalProteins * res.proteinsPassed) / totalDBProteins // should be calculated for each speciesQtyProtein as expected number of
	// proteins attributed to the speciesQtyProtein
	var expNon = (totalProteins * (totalProteins - res.proteinsPassed)) / totalDBProteins
	var pValue = countFisher(totalProteins, m, expSpec, expNon) // Fisher test should be run for each possible speciesQtyProtein,
	// instead of totalProteins should be proteins in the assembly, attributed to the speciesQtyProtein
	fmt.Println(strconv.FormatFloat(pValue, 'E', -1, 32))
	//fmt.Println(time.Now().Format("2006-01-02 15:04:05"))
	//fmt.Println("Finish")
}

func countFisher(totalProteinsFromFile int, specieProteinsFromFile int, expSpec int, expNon int) float64 {
	var pValue, _, _, _ = fisher.FisherExactTest(totalProteinsFromFile, specieProteinsFromFile, expSpec, expNon)
	return pValue
}

func containsString(s []string, searchterm string) bool {
	i := sort.SearchStrings(s, searchterm)
	return i < len(s) && s[i] == searchterm
}

func containsSpecie(species []SpeciesCount, specie string) bool {
	for _, s := range species {
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

func searchProtein(db *sql.DB, hashes <-chan *SeqHash) ([]int, []SpeciesCount) {
	var res []int
	var searchingProteinString string
	for hash := range hashes {
		searchingProteinString = searchingProteinString + "'" + hash.hash + "'" + ", "
	}
	searchingProteinString = strings.TrimSpace(searchingProteinString)
	searchingProteinString = strings.TrimSuffix(searchingProteinString, ",")
	sqlStatement := "select phd.id " +
		"from ProteinHashed as phd " +
		"where phd.Hash in (" + searchingProteinString + ")" +
		"group by phd.id" // this's a way for SQL injections?
	//rows, err := db.Query(sqlStatement, searchingProteinString)
	stmt, err := db.Prepare(sqlStatement)
	rows, err := stmt.Query()
	defer rows.Close()
	if err != nil {
		log.Fatalln(err)
	}
	for rows.Next() {
		var ipgid int
		err = rows.Scan(&ipgid)
		if err != nil {
			log.Fatalln(err)
		}
		res = append(res, ipgid) // should be returned through channel
	}
	var speciesCount []SpeciesCount
	sqlStatement = "SELECT ss.Name,count(phd.id) " +
		"FROM ProteinHashed as phd " +
		"inner join ProteinHashedSpecies as phs on phd.id = phs.ProteinHashed " +
		"inner join Species as ss on phs.Species = ss.id " +
		"where phd.Hash in (" + searchingProteinString + ") " +
		"GROUP by ss.Name"
	stmt, err = db.Prepare(sqlStatement)
	checkError(err)
	defer stmt.Close()
	rows, err = stmt.Query()
	for rows.Next() {
		var specieName string
		var proteinQty int
		err = rows.Scan(&specieName, &proteinQty)
		checkError(err)
		var s SpeciesCount
		s.name = specieName
		s.qty = proteinQty
		speciesCount = append(speciesCount, s)
	}
	return res, speciesCount
}

func searchSpeciesWithProtein(db *sql.DB, protein string) []string {
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
			isDot := id[len(id)-2 : len(id)-1]
			if isDot == "." {
				id = id[:len(id)-2]
			}

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

func writeSpecies(db *sql.DB, species_records map[string][]string) {
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

func countDBProteins(db *sql.DB) int {
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
	return proteinsCount
}
