package main

import (
	"database/sql"
	"fmt"
	"io"
	"os"

	_ "github.com/mattn/go-sqlite3"
	"github.com/shenwei356/bio/seqio/fastx"
)

func main()  {
	file := "c:\\temp\\PAO1.fasta"
	//outfh, err := xopen.Wopen(file)
	//checkError(err)
	//defer outfh.Close()

	reader, err := fastx.NewDefaultReader(file)
	checkError(err)
	//var record *fastx.Record
	for {
		record, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			checkError(err)
			break
		}
		record.Format(0)
		// fmt is slow for output, because it's not buffered
		//fmt.Printf("%s", record.Format(0))
		//record.FormatToWriter(outfh, 0)
	}
	createDB()
}

func checkError(err error) {
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func createDB(){
	if _, err := os.Stat("data.db"); os.IsNotExist(err){
		os.Create("data.db")
	}

	db, err := sql.Open("sqlite3", "data.db")
	checkError(err)
	var tableExists int
	row := db.QueryRow("select count(*) as tableExists from sqlite_master where type='table' AND name='ProteinHashed'")
	row.Scan(&tableExists)
	if tableExists == 0 {
		_, err = db.Exec("CREATE TABLE 'ProteinHashed' ('id' INTEGER PRIMARY KEY, 'Hash' CHAR (100))")
		checkError(err)
	}
	db.Close()
}