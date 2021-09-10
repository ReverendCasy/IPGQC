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
	"strconv"
	"strings"

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

type Records []Record

//TODO: add cmd parameters
var db = flag.String("dbfile", "", "path to database file")

func main() {
	flag.Parse()
	dbName := *db
	db := InitDB(dbName)
	defer db.Close()
	createTable(db)
	stats_file := "test_stats.tsv"
	stats_data := parseTSV(stats_file)
	file := "test_seqs.fasta"
	//outfh, err := xopen.Wopen(file)
	//checkError(err)
	//defer outfh.Close()
	hashes, err := readFastaFile(file)
	if err != nil {
		log.Println(err)
	}
	records := prepareRecords(stats_data, hashes)
	writeToDB(db, records)
	fmt.Println("Finish")
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

func createTable(db *sql.DB) {
	var tableExists int
	row := db.QueryRow("select count(*) as tableExists from sqlite_master where type='table' AND name='ProteinHashed'")
	row.Scan(&tableExists)
	if tableExists == 0 {
		_, err := db.Exec("CREATE TABLE 'ProteinHashed' ('id' INTEGER PRIMARY KEY, 'Hash' CHAR (100))")
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
		checkError(err)
	}
	tx.Commit()
}

func parseTSV(fileName string) map[string][]string {
	tsv, err := os.Open(fileName)
	if err != nil {
		fmt.Println(err)
		return nil
	}

	scn := bufio.NewScanner(tsv)

	var lines []string

	for scn.Scan() {
		line := scn.Text()
		lines = append(lines, line)
	}

	if err := scn.Err(); err != nil {
		fmt.Println(err)
		return nil
	}

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
