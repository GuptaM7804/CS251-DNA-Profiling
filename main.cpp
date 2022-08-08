/*main.cpp*/

//
// Author: Manav Gupta
//
// Project 1: DNA Profilling
//
// Date: 1/25/2022
//
// Creative component: Add user function
// This function can only be called after loading a database as
// this function requires DNA sequences/STR from a database.
// The function is called with inputting "add_user"
// then it requires "name" of an individual
// then it requires numbers for STR repeats
// in the form 2,6,7,...
// The function adds a user of inputted name
// with inputted numbers of DNA sequence repeats.


#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include "ourvector.h"

using namespace std;


// profiler struct- each struct is an element of profilers vector
struct Profilers {
  string name;
  ourvector<int> DNAnum;
};


// Function prototypes
void loadData(string filename, ourvector<Profilers> &profilers, int &x);
void LoadDNAsequence(string &line, ourvector<ourvector<char>> &sequences);
string getName(string &line);
int getDNAnum(string &line);
void loadDataV2(string filename, ourvector<char> &actualDNA);
void processingCounts(ourvector<char> actualDNA,
ourvector<ourvector<char>> sequences, ourvector<int> &numRepeats);
bool matching(ourvector<int> DNAnum, ourvector<int> numRepeats);
void displayProcessedData(ourvector<ourvector<char>> sequences,
ourvector<int> numRepeats);
void load_db(string filename, ourvector<Profilers> &profilers,
int &x, ourvector<ourvector<char>> &sequences);
void IFload_db(ourvector<Profilers> profilers);
void IFload_dna(ourvector<char> actualDNA);
void display(ourvector<Profilers> profilers, int x, int y, int z,
ourvector<ourvector<char>> sequences, ourvector<char> actualDNA,
ourvector<int> numRepeats);
void load_dna(string &filename, ourvector<char> &actualDNA, int &y);
void process(int x, int y, int &z, ourvector<char> actualDNA,
ourvector<int> &numRepeats, ourvector<ourvector<char>> sequences);
void search(int x, int y, int z, ourvector<Profilers> profilers,
ourvector<int> numRepeats);
void add_user(int x, ourvector<Profilers> &profilers,
ourvector<ourvector<char>> sequences);


// This function reads in the data stored in filename and
// fills the ourvector of Profiler structs with that data.
// A Profiler struct contains the name and
// an integer vector of DNA sequence numbers of each profiler.
// It also calls the function that stores give DNA sequences
// into a vector of character vectors-sequences.
void loadData(string filename, ourvector<Profilers> &profilers,
ourvector<ourvector<char>> &sequences, int &x) {
  ifstream inFile(filename);
  string line;

  getline(inFile, line, '\n');


  // remove word "name" from database as it is not necessary
  string tmp = getName(line);
  while (!line.empty()) {
    // loads DNA in vector of vector of char until the line is empty
    LoadDNAsequence(line, sequences);
  }


  // loops until end of file
  while (!inFile.eof()) {
    getline(inFile, line, '\n');
    if (!inFile.fail()) {
      Profilers newProfiler;
      newProfiler.name = getName(line);
      while (!line.empty()) {
        newProfiler.DNAnum.push_back(getDNAnum(line));
      }
      profilers.push_back(newProfiler);
    }
  }
}


// This function reads in and stores the whole DNA sequence
// from a DNA sequence file/DNA to be matched file.
// The sequence is stored in a character vector called actualDNA
void loadDataV2(string filename, ourvector<char> &actualDNA) {
  ifstream inFile(filename);
  string line;

  getline(inFile, line, '\n');
  for (int x = 0; x < line.size(); x++) {
      actualDNA.push_back(line[x]);
    }
}


// This function loads a vector of vectors of characters
// with DNA sequence characters
void LoadDNAsequence(string &line, ourvector<ourvector<char>> &sequences) {
  ourvector<char> temp;


  // if "," is found...
  if (line.find(",") != string::npos) {
    size_t pos = line.find(",");
    for (int x = 0; x < pos; x++) {
      temp.push_back(line[x]);
    }


    // delete elements in line that have already been extracted
    line = line.substr(pos + 1, line.size());
  } else {
    for (int x = 0; x < line.size(); x++) {
      temp.push_back(line[x]);
    }


    // empty line after all elements have been extracted
    line.clear();
  }


  // put each vector of char in a vector to store it all
  sequences.push_back(temp);
}


// This function extracts the name out of the string line and returns it.
string getName(string &line) {
  size_t pos = line.find(",");
  string name = line.substr(0, pos);


  // updates line so the name and "," have been removed
  line = line.substr(pos + 1, line.size() - name.size() - 1);
  return name;
}


// This function extracts all DNA nums listed in the string line and returns it.
int getDNAnum(string &line) {
  size_t pos = line.find(",");
  if (line.find(",") != string::npos) {
    string semStr = line.substr(0, pos);
    line = line.substr(pos+1, line.size() - semStr.size() - 1);
    int strDNA = stoi(semStr);
    return strDNA;
  }
  string semStr = line.substr(0, line.size());
  line.clear();
  int strDNA = stoi(semStr);
  return strDNA;
}


// this function processes the vector of the actual DNA strand given and
// compares it with the DNA sequences in from the database.
// It stores the amount of repeats of each sequence in the actual strand,
// in another vector, of ints, called numRepeats
void processingCounts(ourvector<char> actualDNA,
ourvector<ourvector<char>> sequences, ourvector<int> &numRepeats) {
  for (int i = 0; i < sequences.size(); i++) {
    int count = 0;
    for (int j = 0; j < actualDNA.size() - sequences[i].size(); j++) {
      int ctr = 0;


      // if same character appears in actualDNA as the given sequences...
      if (actualDNA[j] == sequences[i][ctr]) {
        for (int k = j; k < sequences[i].size() + j; k++) {
          // if all the letter match in actualDNA with sequences given...
          if (actualDNA[k] == sequences[i][ctr]) {
            // counter increases
            ctr++;
          }
        }


        // if counter = the size of sequences
        // (ie: all characters matched with given sequences)...
        if (ctr == sequences[i].size()) {
          // count increases
          count++;
        }
      }
    }


    // the count (number of times sequence appears in DNA strand)
    // pushed back into a vector of ints
    numRepeats.push_back(count);
  }
}


// This function returns true if numbers in DNAnum
// match the numbers in numRepeats
bool matching(ourvector<int> DNAnum, ourvector<int> numRepeats) {
  for (int i = 0; i < DNAnum.size(); i++) {
    if (numRepeats[i] != DNAnum[i]) {
      // false returned if they do not match
      return false;
    }
  }


  // the person matched with the given DNA strand
  return true;
}


// this function ouputs the processed data,
// the sequences from the database and
// how many times they appear in DNA strand
void displayProcessedData(ourvector<ourvector<char>> sequences,
ourvector<int> numRepeats) {
  cout << "DNA processed, STR counts: " << endl;
  for (int i = 0; i < sequences.size(); i++) {
    for (int j = 0; j < sequences[i].size(); j++) {
      cout << sequences[i][j];
    }
    cout << ":";
    cout <<  " " << numRepeats[i];
    cout << endl;
  }
  cout << endl;
}


// function that loads the database into profilers and sequences and updates x
void load_db(string filename, ourvector<Profilers> &profilers, int &x,
ourvector<ourvector<char>> &sequences) {
  cout << "Loading database..." << endl;
  cin >> filename;
  ifstream myfile;
  myfile.open(filename);

  if (!myfile.fail()) {
    x = 1;
    loadData(filename, profilers, sequences, x);
    myfile.close();
  } else {
    cout << "Error: unable to open '" << filename << "'" << endl;
  }
}


// function that outputs the loaded database
void IFload_db(ourvector<Profilers> profilers) {
  cout << "Database loaded: " << endl;
  for (int i = 0; i < profilers.size(); i++) {
    cout << profilers[i].name;
    for (int j = 0; j < profilers[i].DNAnum.size(); j++) {
      cout << " " << profilers[i].DNAnum[j];
    }
    cout << endl;
  }
  cout << endl;
}


// function that outputs loaded DNA
void IFload_dna(ourvector<char> actualDNA) {
  cout << "DNA loaded: " << endl;
  for (int i = 0; i < actualDNA.size(); i++) {
    cout << actualDNA.at(i);
  }
  cout << endl;
  cout << endl;
}


// function that collects the output functions above and
// checks whether all the necessary functions have been called or not
void display(ourvector<Profilers> profilers, int x, int y, int z,
ourvector<ourvector<char>> sequences, ourvector<char> actualDNA,
ourvector<int> numRepeats) {
  if (x == 1) {
    IFload_db(profilers);
  } else {
    cout << "No database loaded." << endl;
  }
  if (y == 1) {
    IFload_dna(actualDNA);
  } else {
    cout << "No DNA loaded." << endl;
    cout << endl;
  }
  if (z == 1) {
    displayProcessedData(sequences, numRepeats);
  } else {
    cout << "No DNA has been processed." << endl;
  }
}


// function that opens the DNA strand file and
// loads it into actualDNA vector, also updates y to 1 if file opens
void load_dna(string &filename,ourvector<char> &actualDNA, int &y) {
  cout << "Loading DNA..." << endl;
  cin >> filename;
  ifstream myfile;
  myfile.open(filename);

  if (!myfile.fail()) {
    y = 1;
    loadDataV2(filename, actualDNA);
    myfile.close();
  } else {
    cout << "Error: unable to open '" << filename << "'" << endl;
  }
}


// function that compiles all the processing command
// and checks whether previously necessary functions have been called
// then processes the data
void process(int x, int y, int &z, ourvector<char> actualDNA,
ourvector<int> &numRepeats, ourvector<ourvector<char>> sequences) {
  if (x == 1) {
    if (y == 1) {
      z = 1;
      cout << "Processing DNA..." << endl;
      processingCounts(actualDNA, sequences, numRepeats);
    } else {
      cout << "No DNA loaded." << endl;
    }
  } else {
    cout << "No database loaded." << endl;
  }
}


// function that checks whether previously necessary functions have been called
// then outputs whether DNA matches a profiler or not
void search(int x, int y, int z, ourvector<Profilers> profilers,
ourvector<int> numRepeats) {
  if (x == 1) {
    if (y == 1) {
      if (z == 1) {
        string matched = "";
        cout << "Searching database..." << endl;
        for (int i = 0; i < profilers.size(); i++) {
          if (matching(profilers[i].DNAnum, numRepeats) == true) {
            matched = profilers[i].name;
          }
        }
        if (matched.size() > 0) {
          cout << "Found in database!  DNA matches: " << matched << endl;
        } else {
          cout << "Not found in database" << endl;
        }
      } else {
        cout << "No DNA processed." << endl;
      }
    } else {
      cout << "No DNA loaded." << endl;
    }
  } else {
    cout << "No database loaded." << endl;
  }
}


// adds a user into profilers from user input
void add_user(int x, ourvector<Profilers> &profilers,
ourvector<ourvector<char>> sequences) {
  if (x == 1) {
    Profilers newProfiler;
    cout << "Enter name of profiler: ";
    string username;
    cin >> username;
    cout << "How many repeats of the following DNA does the user have? ";
    cout << "(separate each number with a comma and no spaces 1,4,7)" << endl;
    for (int i = 0; i < sequences.size(); i++) {
      for (int j = 0; j < sequences[i].size(); j++) {
        cout << sequences[i][j];
      }
      cout << " ";
    }
    cout << endl;
    string str;
    cin >> str;
    newProfiler.name = username;
    while (!str.empty()) {
      newProfiler.DNAnum.push_back(getDNAnum(str));
    }
    profilers.push_back(newProfiler);
  }
}


// main function calls all the above function and has the menu of the app
int main() {
  cout << "Welcome to the DNA Profiling Application." << endl;


  // vector of the profilers in database
  ourvector<Profilers> profilers;


  // vector of DNA sequences in database
  ourvector<ourvector<char>> sequences;


  // vector of characters contains all the characters
  // that appear in the DNA strand file
  ourvector<char> actualDNA;


  // vector of integers contains the number of times sequences
  // from databse appear in actualDNA
  ourvector<int> numRepeats;


  // declaring variables for true or false check
  // x checks whether load_db has been called
  // y checks whether load_dna has been called
  // z checks whether process has been called
  int x = 0, y = 0, z = 0;


  // command variable
  string cmd;
  cout << "Enter command or # to exit: ";
  cin >> cmd;
  string filename;


  // while loop, keeps taking in user input until "#" is entered
  while (cmd != "#") {
    if (cmd == "load_db") {
      load_db(filename, profilers, x, sequences);
    } else if (cmd == "display") {
      display(profilers, x, y, z, sequences, actualDNA, numRepeats);
    } else if (cmd == "load_dna") {
      load_dna(filename, actualDNA, y);
    } else if (cmd == "process") {
      process(x, y, z, actualDNA, numRepeats, sequences);
    } else if (cmd == "search") {
      search(x, y, z, profilers, numRepeats);
    } else if (cmd == "add_user") {
      add_user(x, profilers, sequences);
    } else {
      cout << "Command not detected! Try again" << endl;
    }
    cout << "Enter command or # to exit: ";


    // keeps the loop running
    cin >> cmd;
  }
  return 0;
}
