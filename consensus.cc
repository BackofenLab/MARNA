#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stack>
#include <stdlib.h>
#include <string>
#include <vector>

const int DEBUG=0;
const int MIN_BP_DIST=2;
const double CONSERVATION_FACTOR=0.5;
// CONSERVED_CHAR_FACTOR must be greater than 50%
const double CONSERVED_CHAR_FACTOR=0.5;

const double CONSERVED_SINGLESEQ_BASEPAIR_THRESHOLD=0.1;

typedef struct pair {int p1; int p2;} Pair;

using namespace std;

inline void error(const char* msg) {cerr<<"Error: " << msg <<endl;exit(0);}

inline Pair pairpop(stack<Pair> s) { Pair a=s.top(); s.pop(); return a;}
inline int intpop(stack<int> s) { int a=s.top(); s.pop(); return a;}
inline char charpop(stack<char> s) { char a=s.top(); s.pop(); return a;}

// trace entry. we have the following types for an tr(c,c1):
// tr(c,c1).split=false: 
//     tr(c,c1).k=0    =>  consensus_M(c,c1)=consensus_M(c+1,c1-1) + basebp
//     tr(c,c1).k=+1   =>  consensus_M(c,c1)=consensus_M(c+1,c1)
//     tr(c,c1).k=+1   =>  consensus_M(c,c1)=consensus_M(c,c1-1)
// tr(c,c1).split=true =>  consensus_M(c,c1)=consensus_M(c,c+tr(c,c1).k)
//                                          +consensus_M(c+tr(c,c1).k+1,c1)
typedef struct tre {int k; bool split;} TraceEntry;


// ==============================================================//
//                 CLASS Str                                     //
// ==============================================================//
//
class Str {
 protected:

 public:
  string str_string;
  vector<int> str_pos;
  Str();
  Str(string  str_string_);

  friend istream & operator >> (istream &in, Str &str);
  friend ostream & operator << (ostream & out, Str &str);
  void print_str_pos(ostream & out);
  void recalc_bondings();
};




// ==============================================================//
//                 CLASS SeqStr                                  //
// ==============================================================//
//
class SeqStr {
 protected:

 public:
  int num_of_strs;
  SeqStr();
  string seq;
  vector<Str*> strs;
  vector<double> prob;
  void search_for_identifier(istream &in);
  double combined_basepair_prob(int i, int j);
  friend istream & operator >> (istream &in, SeqStr &seqstr);
  friend ostream & operator << (ostream & out, SeqStr &seqstr);

};


// ==============================================================//
//                 CLASS ConsensusStr                            //
// ==============================================================//
//

class ConsensusStr {
 protected:
  bool conserved_cols_calced;
  bool conserved_bp_calced;

 public:
  int num_of_seqs;
  int num_of_cols;
  ConsensusStr();
  vector< vector<int> > alignment; // sorted by rows
  vector<SeqStr*> seqstrs;
  vector<bool> conserved_cols;
  vector<char> conserved_cols_val;
  vector< vector<double> > conserved_basepairs_prob; 

  friend ostream & operator << (ostream & out, ConsensusStr &constr);

  void read_alignment(istream &in);
  void print_alignment_with_strs(ostream &out);
  void print_consensus_seq(ostream &out);
  void read_seqstrs(istream &in);
  void calc_conserved_cols();
  void calc_conserved_basepairs();
  void calc_and_print_consensus_struct(double threshold,ostream &out);
  void backtrace_consensus(int i, int j, 
                           vector< vector<double> > &consensus_M,
                           vector< vector<TraceEntry> > &trace_M,
                           ostream &out);
  char col_print_dot_or_gap(int c);
  
};



// ==============================================================//
//                 CLASS Str                                     //
// ==============================================================//
//
// -------------------------
// Constructor/Destructor
//
Str::Str() {
};

Str::Str(string  str_string_) : str_string(str_string_) {
  recalc_bondings();
};


// -------------------------
// Input/Output
//
istream & operator >> (istream &in, Str &str) {
  in >> str.str_string;
  str.recalc_bondings();
  return in;
}

ostream & operator << (ostream & out, Str &str) {
  out << str.str_string;
 return out;
}

void Str::print_str_pos(ostream & out) {
  for (int i=0;i<(int) str_pos.size();i++)
    out << str_pos[i] << " " ;
}


// -------------------------
// calc-bondings
//
// generates an array of structure positions, stored in str_pos.
// foreach i, str_pos[i]=j if (i,j) are bonded, and -1 else.
//
void Str::recalc_bondings() {

  stack<int> posstack;

  str_pos.resize(str_string.size());

  for (int i=0;i<(int) str_string.size();i++)
    if(str_string[i]=='(')
      posstack.push(i);
    else if(str_string[i]==')')
      {
        if (posstack.empty()) 
          error("unnested structure found ");
        
        str_pos[i]=posstack.top();
        str_pos[posstack.top()]=i;
        posstack.pop();
      }
    else
      str_pos[i]=-1;

  if (!(posstack.empty()))
    error("unnested structure found ");
}



// ==============================================================//
//                 CLASS SeqStr                                  //
// ==============================================================//
//
// -------------------------
// Constructor/Destructor
//
SeqStr::SeqStr() {
  num_of_strs=0;
};





// -------------------------
// Input/Output
//
//

void SeqStr::search_for_identifier(istream &in) {

  string aux="";

  // search for indentifier
  //
  while (aux.find_first_of(">") == string::npos) {
    in >> aux;
  }

}

// input sequence-structure in fasta format
// assumes that either the identifyer, or the sequence is the actual line
// reads until next identifier, or until eof
//
istream & operator >> (istream &in, SeqStr &seqstr) {
  bool calc_prob=false;

  string aux;
  Str* auxstr;

  if (DEBUG) cerr << "========== Read New Sequence-Struct =============:" << endl;
  // read sequences
  //
  in >> aux;
  if (DEBUG) cerr << "First line:" << aux << endl;

  // search for a line without >
  //
  while (aux.find_first_of(">") != string::npos) {
    in >> aux;
    if (DEBUG) cerr << "Search for non->:" << aux << endl;
  }
  seqstr.seq=aux;
  if (DEBUG) cerr << "----> interpreted as sequence" << endl;

  // read structures
  //
  in >> aux;
  if (DEBUG) cerr << "Should be the first structure:" << aux << endl;

  // search for structures. We search for structures (= strings of .() )
  // and a possible following probability. the current aux is always
  // a structure string, or the beginning of the next entry.
  //
  // old: (aux.find_first_of(">") == string::npos)
  while ((aux.find_first_not_of("().") == string::npos) && (in.peek() != EOF))
    {

    // read structure (without end_of_line)
    //
    if (aux.find_first_not_of(".()") == string::npos) {

      // make the structure
      //
      auxstr = new Str(aux);
      seqstr.strs.push_back(auxstr);
      if (DEBUG) cerr << "----> interpreted as struct" << endl;
      
      // read probability or next structure
      //
      in >> aux;
      if (DEBUG) cerr << "next poss struct or prob:" << aux << endl;
      // check, whether this is a probability
      //
      if (aux.find_first_not_of("0123456789.") == string::npos) {
        seqstr.prob.push_back(atof(aux.c_str()));
        if (DEBUG) 
          cerr << "----> interpreted as probability" << endl;
        // we must read the next structure
        in >> aux;
        if (DEBUG) cerr << "next poss struct:" << aux << endl;
      }
      else {
        if (DEBUG) cerr << "----> NO probability" << endl;
        calc_prob=true;
      }
    }
    else {
      if (DEBUG) cerr << "----> interpreted as no struct" << endl;
      in >> aux;
      if (DEBUG) cerr << "next poss struct or seq:" << aux << endl;
    }
  }
  if (DEBUG) cerr << "----> interpreted as end or begin of next seq-struct" << endl;

  seqstr.num_of_strs=seqstr.strs.size();

  if (calc_prob) {
  // calc the probabilities
  //
  seqstr.prob.resize(seqstr.num_of_strs);
  for (int i=0; i < seqstr.num_of_strs; i++) 
    seqstr.prob[i]=1/(float) seqstr.num_of_strs;  
  }

  

  return in;
}

// outputs the sequence and the structures
//
ostream & operator << (ostream & out, SeqStr &seqstr) {

  string aux;
  aux = "Strs: ";

  out << "Seq:  " << seqstr.seq << endl;
  
/*   int prev_prec= out.precision(); */
/*   out.precision(5); */

  for (int i=0; i< (int) seqstr.strs.size(); i++)
    {
      out << endl << aux << *(seqstr.strs[i]) << "   " << seqstr.prob[i];
      aux = "      ";
    }

/*   out.precision(prev_prec); */
  
 return out;
}

// -------------------------
// calculate the combined probability of a base-pair
//
// NOTE: We have extended this approach such that
// if a single basepair-prob with probability over a given threshold is
// set, then we count this as one
//
double SeqStr::combined_basepair_prob(int i, int j) {
  double p=0;
  
  for (int str=0; str < num_of_strs; str++) 
    if (strs[str]->str_pos[i]==j)
      p = p + prob[str];

  // NEW: everything over threshold is counted as 1
  // 
  if (p >= CONSERVED_SINGLESEQ_BASEPAIR_THRESHOLD) 
    p=1;

  return p;
}


// ==============================================================//
//                 CLASS ConsensusStr                            //
// ==============================================================//
//
// -------------------------
// Constructor/Destructor
//
ConsensusStr::ConsensusStr() {
  num_of_seqs=0;
  num_of_cols=0;
};

// -------------------------
// Input/Output
//
//
ostream & operator << (ostream & out, ConsensusStr &constr) {
  
  // print seqstrs
  //
  for (int r=0; r < constr.num_of_seqs ; r++) 
    out << *constr.seqstrs[r] << endl;


  // print align
  //
  out << endl << "Alignment:" << endl;
  //
  for (int r=0; r < constr.num_of_seqs ; r++) {
    for (int c=0; c < constr.num_of_cols ; c++) {
/*       cerr << "Row: " << r  */
/*            << " Col: "<< c */
/*            << " Aligval: " << constr.alignment[r][c] */
/*            << " Seq: " ; */
/*       if (constr.alignment[r][c] != -1) */
/*         cerr << constr.seqstrs[r]->seq[constr.alignment[r][c]]; */
/*       else */
/*         cerr << "_"; */
/*       cerr << endl; */
      
      if (constr.alignment[r][c] != -1)
        out << constr.seqstrs[r]->seq[constr.alignment[r][c]];
      else
        out << "_";
/*       cerr << endl; */
    }
    out << endl;
  }

  // print conservation
  //
  if (constr.conserved_cols_calced)
    for (int c=0; c < constr.num_of_cols ; c++) 
      if (constr.conserved_cols[c]) 
      out << "*";
      else
        out << "_";
  out << endl; 


  if (DEBUG>=1) {
    // print conservation matrix
    //
    cerr << endl << "Basepair Conservation:" << endl;
    cerr << "     ";
    for (int c=0; c < constr.num_of_cols ; c++) 
      if (c < 10) 
        cerr << c << "  ";
      else
        cerr << c << " ";
    cerr << endl;

    for (int c=0; c < constr.num_of_cols ; c++) {
      if (c < 10) 
        cerr << c << ":  ";
      else
        cerr << c << ": ";
      for (int c1=0; c1 < constr.num_of_cols ; c1++) 
        if (c1 < c)
          cerr << "   ";
        else if (constr.conserved_basepairs_prob[c][c1] > 9)
          cerr << constr.conserved_basepairs_prob[c][c1] << " ";
        else
          cerr << constr.conserved_basepairs_prob[c][c1] << "  ";
      cerr << endl;
    }
  }

  return out;
}


// ------------------------- 
// print consensus sequence
//
void ConsensusStr::print_consensus_seq(ostream &out) {
  if (conserved_cols_calced)
    for (int c=0; c < num_of_cols ; c++) 
      out << conserved_cols_val[c];
  else 
    error("Conserved cols must be calculated first --- PRINT_CONSENSUS_SEQ");
  
  out << endl; 
}


// -------------------------
// print alignment with the assoicated structures
//
void ConsensusStr::print_alignment_with_strs(ostream &out) {
// print align
  //
  out << endl << "Alignment with Strs:" << endl;
  //
  for (int r=0; r < num_of_seqs ; r++) {
    // print sequences
    for (int c=0; c < num_of_cols ; c++) {
      if (alignment[r][c] != -1) 
        out << seqstrs[r]->seq[alignment[r][c]];
      else
        out << "_";
    }
    out << endl;
    
    // print structures

    for (int s=0; s < seqstrs[r]->num_of_strs ; s++) {
      for (int c=0; c < num_of_cols ; c++) 
        if (alignment[r][c] != -1) 
          out << (seqstrs[r]->strs[s])->str_string[alignment[r][c]];
        else
          out << "_";
      out << "   " << seqstrs[r]->prob[s] << endl;
    }
  }

}

// Positions of aligned sequences are given in a separate file. 
// File format :
// num_of_seqs num_of_cols
// pos_11 pos_12 pos_13 ...
// ...
// pos_n1 pos_n2 posn_n3 ...

void ConsensusStr::read_alignment(istream &in) 
{
  string aux;
  
  in >> num_of_seqs;
  in >> num_of_cols;
  getline(in,aux);  

  alignment.resize(num_of_seqs);

  for (int r=0; r < num_of_seqs ; r++) {
    alignment[r].resize(num_of_cols);
    for (int c=0; c < num_of_cols ; c++) {
        in >> alignment[r][c];
        // in alignment, everything is numbered 1 to N.
        // but in our program, everything is number 0 to N-1
        // Hence
        if (alignment[r][c] != -1)
          alignment[r][c]--;
    }
  }  
}

// Source sequences and structures 

void ConsensusStr::read_seqstrs(istream &in) 
{
  for (int r=0;r<num_of_seqs;r++) 
    {
      seqstrs.push_back(new SeqStr);    
      in >> *seqstrs[r];
    }
}


// -------------------------
// calculation of conserved cols 
//
void ConsensusStr::calc_conserved_cols() {

  // how many columns must be different from gap for conservation
  //
  int num_seqs_conserved= (int) ceil((double) num_of_seqs*CONSERVATION_FACTOR);
  int num_seqs_conserved_char
    = (int) ceil((double) num_of_seqs*CONSERVED_CHAR_FACTOR);

  conserved_cols.resize(num_of_cols);
  conserved_cols_val.resize(num_of_cols);
                             
  for (int c=0; c < num_of_cols; c++) {

    // calc number of entries different from gap
    //
    int num_nongaps=0;
    int num_A=0,num_C=0,num_G=0,num_U=0;

    for (int r=0; r < num_of_seqs; r++) {
      
      if (alignment[r][c] != -1) {
        num_nongaps++;
        switch( seqstrs[r]->seq[alignment[r][c]] ) {
        case 'A':
        case 'a':
          num_A++;
          break;
        case 'C':
        case 'c':
          num_C++;
          break;
        case 'G':
        case 'g':
          num_G++;
          break;
        case 'U':
        case 'u':
        case 'T':
        case 't':
          num_U++;
          break;
        default:
          error("Unknown character in Input sequence");
          break;
        }
      }
    }

    if (DEBUG >= 3) 
      cerr << " numA: " << num_A 
           << " numC: " << num_C
           << " numG: " << num_G
           << " numU: " << num_U
           << endl;
    
    // check conservation
    //
    conserved_cols[c] = (num_nongaps >= num_seqs_conserved);

    // now check for maximum. We assume that CONSERVED_CHAR_FACTOR 
    // is greater than 50% 
    //
    if (num_A >=  num_seqs_conserved_char)
      if (num_A == num_of_seqs)  
        conserved_cols_val[c]='A';
      else
        conserved_cols_val[c]='a';
    else if (num_C >= num_seqs_conserved_char)
      if (num_C == num_of_seqs)  
        conserved_cols_val[c]='C';
      else
        conserved_cols_val[c]='c';
    else if (num_G >= num_seqs_conserved_char)
      if (num_G == num_of_seqs)  
        conserved_cols_val[c]='G';
      else
        conserved_cols_val[c]='g';
    else if (num_U >= num_seqs_conserved_char)
      if (num_U == num_of_seqs)  
        conserved_cols_val[c]='U';
      else
        conserved_cols_val[c]='u';
    else if (conserved_cols[c])
      conserved_cols_val[c]='N';
    else
      conserved_cols_val[c]='_';
  }
  conserved_cols_calced=true;
}




// -------------------------
// calculation of conserved basepairs. The result is a boolean matrix
// of size num_of_cols x num_of_cols. We will use only the upper 
// triangular part. It calculates the expected number of basepairs seen
// in the different sequences. 
//
void ConsensusStr::calc_conserved_basepairs() {

  if (!conserved_cols_calced) 
        error("Conserved cols must be calculated first --- CALC_CONSERVED_BASEPAIRS");
  else {
    
    conserved_basepairs_prob.resize(num_of_cols);
    for (int c=0; c < num_of_cols; c++) {
      conserved_basepairs_prob[c].resize(num_of_cols);

      if (conserved_cols[c])
        // now for all bigger colums (distance at least 1)
        //
        for (int c1=c+1; c1 < num_of_cols; c1++) {
          if ((conserved_cols[c]) && (conserved_cols[c1])) {
          
            double p=0;
          
            // go through all sequences
            for  (int r=0; r < num_of_seqs; r++) {
              // check, whether sequence is aligned, then add expected
              // value for this. Note that in the alignment, the sequence
              // are numbered 1 to N. But in everything else, 
              //
              if ((alignment[r][c] != -1) && (alignment[r][c1] != -1)) 
                p = p+seqstrs[r]->combined_basepair_prob(alignment[r][c],
                                                         alignment[r][c1]);
            }
            conserved_basepairs_prob[c][c1] = p;
          }
          else
            conserved_basepairs_prob[c][c1] = 0;
        }
    }
    
  }

}




// ---------------------------------
// calculate the consensus structure
// it takes all conserved basepairs (where conserved is defined via
// threshold) and searches for the consensus structure maximizing the
// expected occurrences of basepairs.
//
void ConsensusStr::calc_and_print_consensus_struct(double threshold, 
                                                   ostream &out) {

  double threshold_exp = threshold*num_of_seqs;
  if (DEBUG >= 1)
    cerr << "Threshold: " <<  threshold_exp << " expected number of sequence structures" << endl;

/*   double consensus_M[num_of_cols][num_of_cols]; */
/*   int trace_M[num_of_cols][num_of_cols]; */
/*   stack<Pair> backtrace; */

  vector< vector<double> > consensus_M;
  vector< vector<TraceEntry> > trace_M;
  
  consensus_M.resize(num_of_cols);
  trace_M.resize(num_of_cols);
  for (int c=0; c < num_of_cols; c++) {
    consensus_M[c].resize(num_of_cols);
    trace_M[c].resize(num_of_cols);
  }

  int  split_k=-1;

  double act_bp=-1, max_split=-1;

  // int max_entry=-1, maxi=-1, maxj=-1, i=0, j=0;
  // stack<char> rightout;

  // initialisation
  //
  for (int c=0; c < num_of_cols; c++) {
    consensus_M[c][c]=0;
    for (int k=1; k <= MIN_BP_DIST+1;k++)
      if (c < num_of_cols-k) 
        consensus_M[c][c+k]=0;
  }

  // .............................................
  // start with the smallest distance
  //
  for (int d=MIN_BP_DIST; d < num_of_cols; d++) 
    for (int c=0; c < num_of_cols-d; c++) {
      
      // calculate split value
      //
      max_split=-1;
      split_k=-1;
      for (int k=MIN_BP_DIST; k<d-MIN_BP_DIST; k++) {
        // calc split
        act_bp= consensus_M[c][c+k] + consensus_M[c+k+1][c+d];
        // new max?
        if (act_bp > max_split) {
          max_split=act_bp;
          split_k=k;
        }
      }
        

      // basepairs will be only considered if above treshold
      //
      if (conserved_basepairs_prob[c][c+d]>= threshold_exp) {
        //        cerr << "i: " << c << "j: " << c+d << " is bigger" << endl;
        act_bp=conserved_basepairs_prob[c][c+d] + consensus_M[c+1][c+d-1];
      }
      else
        act_bp=-1;

      // recurions
      consensus_M[c][c+d]
        = max(act_bp, max(consensus_M[c+1][c+d],
                          consensus_M[c][c+d-1]));

      // how was this max achieved generated?
      TraceEntry tr; tr.split=false;
      if (act_bp >= consensus_M[c][c+d])
        tr.k=0;
      else if (consensus_M[c+1][c+d] >= consensus_M[c][c+d])
        tr.k=+1;
      else
        tr.k=-1;
      

      // max by split?
      if (max_split >= consensus_M[c][c+d]) {
        consensus_M[c][c+d]=max_split;
        tr.split=true ; 
        tr.k=split_k;
      }
      
      // trace entry
      trace_M[c][c+d]=tr;
    }
  
    


  // .............................................
  // print consensus matrix
  //
  if (DEBUG>=2) {

    int prev_prec= cerr.precision(); 
    cerr.precision(2); 

    cerr << endl << "Consensus Matrix:" << endl;
    cerr << "      ";
    for (int c=0; c < num_of_cols ; c++) 
      if (c < 10) 
        cerr << c << "   ";
      else
        cerr << c << "  ";
    cerr << endl;

    for (int c=0; c < num_of_cols ; c++) {
      if (c < 10) 
        cerr << c << ":  ";
      else
        cerr << c << ": ";
      for (int c1=0; c1 < num_of_cols ; c1++) 
        if (c1 < c)
          cerr << "   ";
        else if (consensus_M[c][c1] > 9)
          cerr << consensus_M[c][c1] << " ";
        else
          cerr << consensus_M[c][c1] << "  ";
      cerr << endl;
    }

    cerr.precision(prev_prec); 

/*     cerr << "Maximal Value: " << max_entry << " at i: " << maxi  */
/*          << " and j: " << maxj << endl; */
  }

  // .............................................
  // backtrace
  //
  //  cerr << " hier" << endl;
  backtrace_consensus(0,num_of_cols -1, consensus_M, trace_M, out);
  out << endl;
  // cerr << " hier-done" << endl;
}

// ---------------------------------
// backtrace for the calculation of consensus structure
//
void ConsensusStr::backtrace_consensus(int i, int j,
                                       vector< vector<double> > &consensus_M,
                                       vector< vector<TraceEntry> > &trace_M,
                                       ostream &out)
{

  if (j-i >= MIN_BP_DIST) {
    
    // cerr << "i: " << i << " j: " << j << " " << "hier-a" << endl;
    
    // maximum by split?
    //
    if (trace_M[i][j].split) {
      // cerr << "i: " << i << " j: " << j << " " << "hier-b" << endl;
      backtrace_consensus(i,i+trace_M[i][j].k, consensus_M, trace_M,out);
      // cerr << "i: " << i << " j: " << j << " " << "hier-b2" << endl;
      backtrace_consensus(i+trace_M[i][j].k+1,j, consensus_M, trace_M,out);
      // cerr << "i: " << i << " j: " << j << " " << "hier-b3" << endl;
    }
    else if (trace_M[i][j].k==0) {
      // maxum by pairing
      // cerr << "i: " << i << " j: " << j << " " << "hier-d" << endl;

        // cerr << "i: " << i << " j: " << j << " " << "hier-d1" << endl;
        out << '('; 
        backtrace_consensus(i+1,j-1, consensus_M, trace_M,out);
        // cerr << "i: " << i << " j: " << j << " " << "hier-d1-done" << endl;
        out << ')';
      }
    else if (trace_M[i][j].k==+1) {
        // cerr << "i: " << i << " j: " << j << " " << "hier-d2" << endl;
        out << col_print_dot_or_gap(i); 
        backtrace_consensus(i+1,j, consensus_M, trace_M,out);
        // cerr << "i: " << i << " j: " << j << " " << "hier-d2-done" << endl;
    }
    else {
        // cerr << "i: " << i << " j: " << j << " " << "hier-d3" << endl;
        backtrace_consensus(i,j-1, consensus_M, trace_M,out);
        // cerr << "i: " << i << " j: " << j << " " << "hier-d3-done" << endl;
        out << col_print_dot_or_gap(j); 
    }
  }
  else {
    // no basepairing possible
    //
    // cerr << "i: " << i << " j: " << j << " " << "hier-e" << endl;
    for (int k=0; k <= j-i; k++) 
      out << col_print_dot_or_gap(i+k);
    //cerr << "i: " << i << " j: " << j << " " << "hier-e done" << endl;
  }

}
    
// -------------------    
// auxilary function for printing consensus_str.
//
char ConsensusStr::col_print_dot_or_gap(int c) {
  if (conserved_cols[c])
    return '.';
  else
    return '_';
}


int main(int argc, char **argv) {

  // Str str1;
  //  SeqStr* seqstr1 = new SeqStr;
  ConsensusStr c;

  char *align_file;
  char *seqstr_file;
  int  gap_length=0;

  if (argc >= 4)
    {
      align_file =argv[1];
      seqstr_file=argv[2];
      gap_length =atoi(argv[3]);
    }
  else {
      align_file=(char*)"positions.aln";
      seqstr_file=(char*)"formatted_sequences.efasta";
  }

  ifstream in_align(align_file);
  ifstream in_seqstr(seqstr_file);

  c.read_alignment(in_align);
  c.read_seqstrs(in_seqstr);
  in_align.close();  
  in_seqstr.close();  
  c.calc_conserved_cols();
  c.calc_conserved_basepairs();

  cout << "Seq      " ;
  for(int i=0;i<gap_length-9;i++) 
    cout << ' ';
  c.print_consensus_seq(cout);

  for(double i=0.3;i<=1.0;i+=0.1) 
    {
      cout << "Str(" << 100*i << "%)" ;
      if (100*i<99) cout << " ";
      for(int j=0;j<gap_length-9;j++) 
	cout << ' ';
      c.calc_and_print_consensus_struct(i,cout);
    }
  return(0);
}
