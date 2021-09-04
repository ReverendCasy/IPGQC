package main

import (
	"bufio"
	"crypto/md5"
	"database/sql"
	"encoding/hex"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	_ "github.com/mattn/go-sqlite3"
	"github.com/shenwei356/bio/seqio/fastx"
)

type  Record struct {
	id int
	hash string
}

type Records []Record

func main()  {
	dbName := "data.db"
	db := InitDB(dbName)
	defer db.Close()
	createTable(db)
	stats_file := "test_stats.tsv"
	stats_data := parseTSV(stats_file)
	file := "test_seqs.fasta"
	//outfh, err := xopen.Wopen(file)
	//checkError(err)
	//defer outfh.Close()
	var r Records
	reader, err := fastx.NewDefaultReader(file)
	checkError(err)
	//var record *fastx.Record
	counter := 0
	for {
		counter += 1
		record, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			checkError(err)
			break
		}
		// fmt is slow for output, because it's not buffered
		//fmt.Printf("%s", record.Format(0))
		id := string(record.ID)
		id = id[:len(id)-2]
		ipgId := stats_data[id]
		hash := getMD5Hash(string(record.Seq.Seq))
		s, err := strconv.Atoi(ipgId[1])
		checkError(err)
		rec := Record{s, hash}
		r = append(r, rec)
		//fmt.Printf("%s", record.Seq.Seq)
		//record.FormatToWriter(outfh, 0)
	}
	writeToDB(db, r)
	fmt.Println("Finish")
}

func checkError(err error) {
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func  InitDB(dbName string) *sql.DB {
	if _, err := os.Stat(dbName); os.IsNotExist(err){
		os.Create(dbName)
	}
	db, err := sql.Open("sqlite3", dbName)
	checkError(err)
	return db
}

func createTable(db *sql.DB){
	var tableExists int
	row := db.QueryRow("select count(*) as tableExists from sqlite_master where type='table' AND name='ProteinHashed'")
	row.Scan(&tableExists)
	if tableExists == 0 {
		_, err := db.Exec("CREATE TABLE 'ProteinHashed' ('id' INTEGER PRIMARY KEY, 'Hash' CHAR (100))")
		checkError(err)
	}
}

func writeToDB(db *sql.DB, items []Record)  {
	sqlStatement := "insert into ProteinHashed (id, Hash) values ($1, $2)"
	stmt, err := db.Prepare(sqlStatement)
	checkError(err)
	defer stmt.Close()
	tx, err := db.Begin()
	checkError(err)
	for _, item := range items{
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