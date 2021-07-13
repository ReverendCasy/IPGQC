//#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <locale>
#include <set>
#include "md5.h"
#include "md5.cpp"
#include <sqlite3.h>


using namespace std;

// Sql callback
static int callback(void *data, int argc, char **argv, char **azColName){
    int i;
    fprintf(stderr, "%s: ", (const char*)data);

//    for(i = 0; i<argc; i++){
//        printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
//    }
    return 0;
}


int main(int argc, char *argv[]) {

    bool to_make_db = false;
    string path_to_db;
    bool to_make_a_search = false;
    string path_to_names = argv[1];
    string path_to_query;

    // Attributes parsing
    for (int i = 2; i < argc; i++) {
        char a[] = "-d";
        char *attribute = a;
        int co = strcmp(argv[i], attribute);
        if (co == 0) {
            to_make_db = true;
            path_to_db = argv[i + 1];
        }
        char a2[] = "-s";
        char *attribute2 = a2;
        co = strcmp(argv[i], attribute2);
        if (co == 0) {
            to_make_a_search = true;
            path_to_query = argv[i + 1];
            cout << path_to_query << endl;
        }
    }

    fstream db, names, q;

    // Sql variables
    sqlite3 *sql_db;
    char *zErrMsg = 0;
    int rc;
    char *sql;

    // Open sql database
    rc = sqlite3_open("hash_database.db", &sql_db);
    if (rc) {
        fprintf(stderr, "Can't open sql database: %s\n", sqlite3_errmsg(sql_db));
        return(0);
    }
    else {
        fprintf(stdout, "Sql database opened successfully\n");
    }

    // Open files
//    db.open("/home/ruslan/PycharmProjects/Hash_test/test_data/smu_db.fasta"); // simple db
//    db.open("/home/ruslan/PycharmProjects/Hash_test/test_data/all_prokaryotes_ipg.fasta"); // big db
//    names.open("/home/ruslan/PycharmProjects/Hash_test/test_data/smu_ids.tsv"); // table id - protein
//    q.open("/home/ruslan/PycharmProjects/Hash_test/test_data/GCA_000829155.1.fasta"); // выборка


    names.open(path_to_names); // table id - protein

    // Success if everything is open
//    if ((db) && (names) && (q)) {
//        fprintf(stdout, "Files opened successfully\n");
//    }
//    else {
//        fprintf(stderr, "Can't one or more files: %s\n", sqlite3_errmsg(sql_db));
//        return(0);
//    }

    // Read file "names" (id - protein)
    string line, field;
    vector< vector<string> > id_protein;  // the 2D array
    vector<string> v;                // array of values for one line only
    while (getline(names, line)) {
        v.clear();
        stringstream ss(line);
        while (getline(ss,field,'\t')) { // break line into comma delimitted fields
            v.push_back(field);  // add each field to the 1D array
        }
        id_protein.push_back(v);  // add the 1D array to the 2D array
    }

    // Now we have an array id - protein

    //    Print an array id - protein
//    for (size_t i=0; i<id_protein.size(); ++i) {
//        for (size_t j=0; j<id_protein[i].size(); ++j) {
//            cout << id_protein[i][j] << ' ';
//        }
//        cout << "\n";
//    }


    // Option 1: Making and filling a database
    if (to_make_db) {
        db.open(path_to_db); // simple db
        // Create a table in sql database
        sql = "CREATE TABLE HASHED("  \
      "ID INT PRIMARY KEY     NOT NULL," \
      "HASH        CHAR(50));";

        rc = sqlite3_exec(sql_db, sql, callback, 0, &zErrMsg);
        if (rc != SQLITE_OK) {
            fprintf(stderr, "SQL error: %s\n", zErrMsg);
            sqlite3_free(zErrMsg);
            return(0);
        }
        else {
            fprintf(stdout, "Table created successfully\n");
        }


        // Read protein database file
        string name, content;
        int length_of_db = 0;
        while (getline(db, line).good()) {
            v.clear();
            if (line.empty() || line[0] == '>') { // Identifier marker
                if (!name.empty()) { // Print out what we read from the last entry
                    size_t dot_pos = name.find(".");
                    string cut_name = name.substr(0, dot_pos);
                    string hashed_content = md5(content);
                    string id = "NOT FOUND";
                    int flag = 0;
                    for (size_t p = 0; p < id_protein.size(), flag == 0; ++p) {
                        if (cut_name == id_protein[p][1]) {
                            id = id_protein[p][0];
                            flag = 1;
                        }
                    }

                    if (id != "NOT FOUND") {
                        string s = "INSERT INTO HASHED (ID,HASH) VALUES (" + id +  ", '" + hashed_content + "'); ";
                        sql = &s[0];
                        rc = sqlite3_exec(sql_db, sql, callback, 0, &zErrMsg);
                        length_of_db++;
                    }
                    name.clear();
                }
                if (!line.empty()) {
                    name = line.substr(1);
                }
                content.clear();
            } else if (!name.empty()) {
                if (line.find(' ') != string::npos) { // Invalid sequence--no spaces allowed
                    name.clear();
                    content.clear();
                } else {
                    content += line;
                }
            }
        }
        if (!name.empty()) { // Print out what we read from the last entry

            size_t dot_pos = name.find(".");
            string cut_name = name.substr(0, dot_pos);
            string hashed_content = md5(content);
            string id = "NOT FOUND";
            int flag = 0;
            for (size_t p = 0; p < id_protein.size(), flag == 0; ++p) {
                if (cut_name == id_protein[p][1]) {
                    id = id_protein[p][0];
                    flag = 1;
                }
            }

            if (id != "NOT FOUND") {
                string s = "INSERT INTO HASHED (ID,HASH) VALUES (" + id +  ", " + hashed_content + "); ";
                sql = &s[0];
                rc = sqlite3_exec(sql_db, sql, callback, 0, &zErrMsg);
                length_of_db++;
            }
        }

        cout << length_of_db << " proteins are successfully added to sql db\n";

        // Now we have an sql db
        db.close();
    }


    // Print the database
//    const char* data = "Callback function called";
//    sql = "SELECT * from HASHED";
//    rc = sqlite3_exec(sql_db, sql, callback, (void*)data, &zErrMsg);
//
//    if( rc != SQLITE_OK ) {
//        fprintf(stderr, "SQL error: %s\n", zErrMsg);
//        sqlite3_free(zErrMsg);
//    } else {
//        fprintf(stdout, "Operation done successfully\n");
//    }


    // Option 2: Searching
    if (to_make_a_search) {
        q.open(path_to_query); // выборка

        // Open the query file, preparing a set for the query
        set <string> st;
        string name, content;

        while (getline(q, line).good()) {
            v.clear();
            if (line.empty() || line[0] == '>') { // Identifier marker
                if (!name.empty()) { // Print out what we read from the last entry
                    size_t dot_pos = name.find(".");
                    string cut_name = name.substr(0, dot_pos);
                    st.insert(md5(content));
                    name.clear();
                }
                if (!line.empty()) {
                    name = line.substr(1);
                }
                content.clear();
            } else if (!name.empty()) {
                if (line.find(' ') != string::npos) { // Invalid sequence--no spaces allowed
                    name.clear();
                    content.clear();
                } else {
                    content += line;
                }
            }
        }
        if (!name.empty()) { // Print out what we read from the last entry
            size_t dot_pos = name.find(".");
            string cut_name = name.substr(0, dot_pos);
            st.insert(md5(content));
        }

        set <string> st_found;


        int counter = 0;
        int total = 0;

        set <string> :: iterator it = st.begin();

        for (int i = 1; it != st.end(); i++, it++) {
            string this_hash = *it;
            string s = "SELECT EXISTS (SELECT 1 FROM HASHED WHERE HASH = '" + this_hash + "')";
            sql = &s[0];
            rc = sqlite3_exec(sql_db, sql, 0, 0, &zErrMsg);
            if (rc == false) {
//                sqlite3_stmt* stmt;
//                stringstream ss;
//                ss << "select * from hashed where hash = '" << this_hash << "';";
//                string sqll(ss.str());
//                sqlite3_prepare_v2(sql_db, sqll.c_str(), -1, &stmt, NULL);
//                sqlite3_finalize(stmt);
//                int ret_code = 0;
//                while((ret_code = sqlite3_step(stmt)) == SQLITE_ROW) {
//                    st_found.insert(sqlite3_column_double(stmt, 0));
//                }

                st_found.insert(this_hash);
                counter++;
            }

            total = i;
        }

        set <string> :: iterator it_found = st_found.begin();

        cout << counter << " of " << total << " found in the database\n";
        cout << "The present IPGs are:\n";

        for (int i = 1; it_found != st_found.end(); i++, it_found++) {
            cout << *it_found << "\n";
        }

        q.close();
        sqlite3_close(sql_db);
    }


    return 1;
}