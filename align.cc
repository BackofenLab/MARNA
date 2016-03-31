#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <list>
#include <math.h>
#include <stack>
#include <string>
#include <vector>
#include <getopt.h>

using namespace std;

string               seq1,seq2,str1,str2;     // original sequences
string               aseq1,aseq2,astr1,astr2; // aligned sequences
vector<double>       weights;
vector<int>          L1,R1,I1,L2,R2,I2;       // L,R:left,right arc positions, I: numbering of arcs
stack<int>           str_stack;
vector<vector<double> > M;
vector<vector<double> > D;
const double eps=0.0000001;
bool tcoffee=false;

// four-tuples are two tuples defining sequence ranges l1,...,r1
// and l2,...,r2 in sequence one and two, respectively.
struct four_tuple{
  int l1;
  int r1;
  int l2;
  int r2;
};

struct alignment{
  int p1, p2;
  double weight;
  alignment* next;
};

double w_d =2.0;  // base deletion
double w_m =1.0;  // base mismatch
double w_r =2.0;  // arc  removing
double w_b =1.5;  // arc  breaking
double w_am=1.8;  // arc  mismatch

int not_free1    (int pos)          { return (str1[pos]=='.' ? 0:1)      ; }
int not_free2    (int pos)          { return (str2[pos]=='.' ? 0:1)      ; }
int base_mismatch(int pos1,int pos2){ return (seq1[pos1]!=seq2[pos2]?1:0); }

double min4(double a,double b,double c,double d)
{
  if (a<=b && a<=c && a<=d)
    return a;
  else if (b<=a && b<=c && b<=d)
    return b;
  else if (c<=a && c<=b && c<=d)
    return c;
  return d;
}

void matrix_output()
{
  cout << "D: " << endl;
  for(int i=0;i<D.size();i++)
    {
      for(int j=0;j<D[i].size();j++)
	cout << D[i][j] << " ";
      cout << endl;
    }
  
  cout << "M: " << endl;
  for(int i=0;i<M.size();i++)
    {
      for(int j=0;j<M[i].size();j++)
	cout << M[i][j] << " ";
      cout << endl;
    }
}

void insert(alignment* start,int p1,int p2,double weight)
{
  // Construction of the alignment
  alignment* insertion = new alignment;
  insertion->p1     = p1;
  insertion->p2     = p2;
  insertion->weight = weight;
  insertion->next   = NULL;
  
  alignment* iter   = start;
  
  if (p1==-1 && p2!=-1)
    while (iter->next!=NULL && p2>iter->next->p2)
      iter=iter->next;
  else if (p2==-1 && p1!=-1)
    while (iter->next!=NULL && p1>iter->next->p1)
      iter=iter->next;
  else
    while (iter->next!=NULL && p2>iter->next->p2 && p1>iter->next->p1 )
      iter=iter->next;
  
  insertion->next=iter->next;
  iter->next=insertion;
}

  /*
    
1st structure :   (         (                  )         )
                 L1[i]  L1[I1[a1]]          L1[i]+k    R1[i]
                           =i1              =a1,k>0

2nd structure :   (         (                  )         )
                 L2[j]  L2[I2[a2]]          L2[j]+l    R2[j]
                           =i2              =a2,l>0


used matrices : 
                 M : "sequence alignment" matrix from L1[i]+1 to R1[i]-1
		 D : base pair alignment matrix, 
		     D(m1,m2) is alignment score from L1[i]+1 to R1[i]-1 
		                                  and L2[j]+1 to R2[j]-1 
		     base pair score must be added outside.
   
  */

double computation()
{
  // computes and returns the alignment score subject to distance values
  D.resize(L1.size());
  for(int s=0;s<D.size();s++)
    D[s].resize(L2.size());
  
  double v1,v2,v3,v4;

  for (int i=0;i<L1.size();i++)
    {
      if(!tcoffee)
	cout << "\r" << 100*i/(L1.size()-1) << "% " << flush ;
      for (int j=0;j<L2.size();j++)
	{
	  // init M
	  M.resize(R1[i]-L1[i]);
	  for(int s=0;s<M.size();s++)
	    M[s].resize(R2[j]-L2[j]);
	  
	  M[0][0]=0;
	  for (int k=1;k<R1[i]-L1[i];k++)
	    M[k][0]=M[k-1][0]+w_d+not_free1(L1[i]+k)*(0.5*w_r-w_d);
	  
	  for (int l=1;l<R2[j]-L2[j];l++)
	    M[0][l]=M[0][l-1]+w_d+not_free2(L2[j]+l)*(0.5*w_r-w_d);
	  
	  // compute M
	  for (int k=1;k<R1[i]-L1[i];k++)
	    for (int l=1;l<R2[j]-L2[j];l++)
	      {
		v1=v2=v3=v4=10000;
		int a1=L1[i]+k;
		int a2=L2[j]+l;
		
		v1=M[k-1][l]+w_d+not_free1(a1)*(0.5*w_r-w_d);
		v2=M[k][l-1]+w_d+not_free2(a2)*(0.5*w_r-w_d);
		v3=M[k-1][l-1]+base_mismatch(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*0.5*w_b;
		
		if (str1[a1]==')' && str2[a2]==')') 
		  {
		    int i1=L1[I1[a1]];
		    int j1=L2[I2[a2]];
		    v4=M[i1-L1[i]-1][j1-L2[j]-1]+D[I1[a1]][I2[a2]]+
		      (base_mismatch(L1[I1[a1]],L2[I2[a2]])+
		       base_mismatch(R1[I1[a1]],R2[I2[a2]]))*0.5*w_am;
		  }
		
		M[k][l]=min4(v1,v2,v3,v4);
	      }
	  D[i][j]=M[R1[i]-L1[i]-1][R2[j]-L2[j]-1];
	}
    }
  
  if(!tcoffee)
    cout << "\r" ;  
  // last run

  // init M
  M.resize(str1.size()+1);
  for(int s=0;s<M.size();s++)
    M[s].resize(str2.size()+1);
  
  M[0][0]=0;
  for (int k=1;k<=str1.size();k++)
    M[k][0]=M[k-1][0]+w_d+not_free1(k-1)*(0.5*w_r-w_d);
  
  for (int l=1;l<=str2.size();l++)
    M[0][l]=M[0][l-1]+w_d+not_free2(l-1)*(0.5*w_r-w_d);
  
  // compute M
  for (int k=1;k<=str1.size();k++)
    for (int l=1;l<=str2.size();l++)
      {
	v1=v2=v3=v4=10000;
	v1=M[k-1][l]+w_d+not_free1(k-1)*(0.5*w_r-w_d);
	v2=M[k][l-1]+w_d+not_free2(l-1)*(0.5*w_r-w_d);
	v3=M[k-1][l-1]+base_mismatch(k-1,l-1)*w_m+(not_free1(k-1)+not_free2(l-1))*0.5*w_b;
	
	if (str1[k-1]==')' && str2[l-1]==')') 
	  {
	    int i1=L1[I1[k-1]];
	    int j1=L2[I2[l-1]];
	    v4=M[i1][j1]+D[I1[k-1]][I2[l-1]]+
	      (base_mismatch(L1[I1[k-1]],L2[I2[l-1]])+
	       base_mismatch(R1[I1[k-1]],R2[I2[l-1]]))*0.5*w_am;
	  }
	M[k][l]=min4(v1,v2,v3,v4);
      }
  
  return M[str1.size()][str2.size()];
}

void traceback()
{
  // stores aligned sequences and weights in aseq1,aseq2,astr1,astr2
  alignment* ali=new alignment;
  ali->p1=-1;  
  ali->p2=-1; 
  ali->next=NULL;
  stack<double> weight;
  double v1,v2,v3,v4;

  // range is the currently computed sequence range
  four_tuple range;
  range.l1=0;
  range.l2=0;
  range.r1=seq1.size()-1;
  range.r2=seq2.size()-1;

  stack<four_tuple> ranges;
  ranges.push(range);
  
  while(!ranges.empty())
    {
      int l1=ranges.top().l1;
      int r1=ranges.top().r1;
      int l2=ranges.top().l2;
      int r2=ranges.top().r2;
      ranges.pop();

      if (l1>r1 && l2<=r2)
	for (int s=r2;s>=l2;s--)
	  insert(ali,-1,s,w_d);
      else if (l1<=r1 && l2>r2)
	for (int s=r1;s>=l1;s--)
	  insert(ali,s,-1,w_d);
      else if (l1<=r1 && l2<=r2)
	{
	  // init and compute M
	  M.resize(r1-l1+2);
	  for(int s=0;s<M.size();s++)
	    {
	      M[s].resize(r2-l2+2);
	      for(int t=0;t<M[s].size();t++)
		M[s][t]=0;
	    }
	  
	  M[0][0]=0;
	  for (int k=1;k<r1-l1+2;k++)
	    M[k][0]=M[k-1][0]+w_d+not_free1(l1+k-1)*(0.5*w_r-w_d);
	  
	  for (int l=1;l<r2-l2+2;l++)
	    M[0][l]=M[0][l-1]+w_d+not_free2(l2+l-1)*(0.5*w_r-w_d);
	  
	  for (int k=1;k<r1-l1+2;k++)
	    for (int l=1;l<r2-l2+2;l++)
	      {
		v1=v2=v3=v4=10000;
		int a1=l1+k-1;                      // a1,a2 sequence positions 
		int a2=l2+l-1;
		
		v1=M[k-1][l]+w_d+not_free1(a1)*(0.5*w_r-w_d);
		v2=M[k][l-1]+w_d+not_free2(a2)*(0.5*w_r-w_d);
		v3=M[k-1][l-1]+base_mismatch(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*0.5*w_b;
		
		if (str1[a1]==')' && str2[a2]==')') 
		  {
		    int i1=L1[I1[a1]];
		    int j1=L2[I2[a2]];
		    v4=M[i1-l1][j1-l2]+D[I1[a1]][I2[a2]]+
		      (base_mismatch(i1,j1)+base_mismatch(a1,a2))*0.5*w_am;
		  }
		M[k][l]=min4(v1,v2,v3,v4);
	      }
	  
	  bool seqaln=true;
	  int k=r1-l1+1;
	  int l=r2-l2+1;
	  // sequence alignment
	  while (seqaln)
	    {
	      int a1=l1+k-1;                      // a1,a2 sequence positions 
	      int a2=l2+l-1;

	      if (k==0 && l==0)
		seqaln=false;
	      else if (k>0 && fabs(M[k][l]-(M[k-1][l]+w_d+not_free1(a1)*(0.5*w_r-w_d)))<eps )
		{
		  insert(ali,a1,-1,w_d+not_free1(a1)*(0.5*w_r-w_d));
		  k--;
		}
	      else if (l>0 && fabs(M[k][l]-(M[k][l-1]+w_d+not_free2(a2)*(0.5*w_r-w_d)))<eps )
		{
		  insert(ali,-1,a2,w_d+not_free2(a2)*(0.5*w_r-w_d));
		  l--;
		}
	      else if (k>0 && l>0 && fabs(M[k][l]-(M[k-1][l-1]+base_mismatch(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*0.5*w_b))<eps)
		{
		  insert(ali,a1,a2,base_mismatch(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*0.5*w_b);
		  k--;
		  l--;
		}
	      else
		seqaln=false;
	    }
	  
	  int a1=l1+k-1;                      // a1,a2 sequence positions 
	  int a2=l2+l-1;                      // right arc ends

	  // base-pair alignment
	  if (str1[a1]==')' && str2[a2]==')')
	    {
	      double w=M[L1[I1[a1]]-l1][L2[I2[a2]]-l2]+D[I1[a1]][I2[a2]]+(base_mismatch(L1[I1[a1]],L2[I2[a2]])+base_mismatch(a1,a2))*0.5*w_am;
	      if (fabs(M[k][l]-w)<eps)
		{
		  int i1=L1[I1[a1]];              // left arc ends
		  int j1=L2[I2[a2]];
		  
		  double edge_weight=0.5*(base_mismatch(L1[I1[a1]],L2[I2[a2]])+base_mismatch(a1,a2))*0.5*w_am;

		  insert(ali,i1,j1,edge_weight);
		  insert(ali,a1,a2,edge_weight);
		  
		  four_tuple CR1,CR2;
		  CR1.l1=l1   ; CR1.r1=i1-1 ; CR1.l2=l2   ; CR1.r2=j1-1 ;
		  CR2.l1=i1+1 ; CR2.r1=a1-1 ; CR2.l2=j1+1 ; CR2.r2=a2-1 ;
		  ranges.push(CR1);
		  ranges.push(CR2);
		}
	    }
	}
    }
  // write aligned sequences
  weights.resize(0);
  for(alignment* iter=ali->next;iter!=NULL;iter=iter->next)
    {
      aseq1.push_back((iter->p1==-1?'-':seq1[iter->p1]));
      astr1.push_back((iter->p1==-1?'-':str1[iter->p1]));
      aseq2.push_back((iter->p2==-1?'-':seq2[iter->p2]));
      astr2.push_back((iter->p2==-1?'-':str2[iter->p2]));
      weights.push_back(iter->weight);
    }
}

void read_parameters(char* filename)
{
  ifstream IFS(filename);
  if (filename!="")
    {
      IFS.seekg(0);
      IFS >> w_d >> w_m >> w_r >> w_b >> w_am;
    }
}

void output_tcoffee()
{
  int p1=1,p2=1;
  for(int i=0;i<aseq1.size();i++)
    {
      if (aseq1[i]!='-')
	p1++;
      if (aseq2[i]!='-')
	p2++;
      // find all realized edges, transform dist into sim according to per cent values
      if (aseq1[i]!='-' && aseq2[i]!='-')
	cout << p1-1 << " " << p2-1 << " " << (int)(100-4*weights[i]+0.5) << endl;
    }
}

void print_usage()
{
  cerr << "\nusage : ./align ACGUACGU \"((....))\" ACGUACGU \".(....).\" " << endl;
  cerr << "Implementation of Jiang's RNA alignment algorithm" << endl;
  cerr << "options"                                                        << endl;
  cerr << "  h      : help"                                                                  << endl;
  cerr << "  t      : t-coffee output, i.e. weighted, realized alignment edges."             << endl;
  cerr << "  s file : read file with scores in this order (default values in parentheses) :" << endl;
  cerr << "            base deletion (w_d =1.0) " << endl;
  cerr << "            base mismatch (w_m =1.0) " << endl;
  cerr << "            arc  removing (w_r =2.0) " << endl;
  cerr << "            arc  breaking (w_b =1.5) " << endl;
  cerr << "            arc  mismatch (w_am=1.8) " << endl; 
}

int main(int argc,char** argv)
{
  int opt;
  while ((opt=getopt(argc,argv,"hts:"))!=-1)
    switch(opt)
      {
      case 'h':
	print_usage();
	exit(0);
      case 't':
	tcoffee=true;
        break;
      case 's':
	read_parameters(optarg);
        break;
      case '?' :
        cerr << "Unrecognized option." << endl;
        exit(1);
        break;
      }

  // Input recognition
  vector<string> input;
  for(;optind<argc;optind++)
    input.push_back(argv[optind]);

  if (input.size()!=4)
    {
      print_usage();
      exit(1);
    }

  // read sequences and structures
  seq1=input[0];
  str1=input[1];
  seq2=input[2];  
  str2=input[3];

  I1.resize(str1.size());
  I2.resize(str2.size());

  // read arcs for sequence 1 
  int index=0;
  for(int i=0;i<str1.size();i++)
    {
      I1[i]=-1;
      if (str1[i]=='(')
	str_stack.push(i);
      else if (str1[i]==')')
	{
	  int last=str_stack.top();
	  str_stack.pop();

	  I1[i]=I1[last]=index++;
	  L1.push_back(last);
	  R1.push_back(i);
	}
    }

  // read arcs for sequence 2
  index=0;
  for(int i=0;i<str2.size();i++)
    {
      I2[i]=-1;
      if (str2[i]=='(')
	str_stack.push(i);
      else if (str2[i]==')')
	{
	  int last=str_stack.top();
	  str_stack.pop();

	  I2[i]=I2[last]=index++;
	  L2.push_back(last);
	  R2.push_back(i);
	}
    }

  if(!tcoffee)
    {
      cout << "used score scheme : (w_d,w_m,w_r,w_b,w_am)=(" ;
      cout << w_d<<','<<w_m<<','<<w_r<<','<<w_b<<','<<w_am<<')'<<endl;
      cout << seq1 << endl;  
      cout << str1 << endl;
      cout << seq2 << endl;  
      cout << str2 << endl;
    }

  double score=computation();
  traceback();

  if(!tcoffee)
    cout << aseq1 << "\n" << astr1 << "\n" << aseq2 << "\n" << astr2 << "\nscore:" << score << endl;
  else
    output_tcoffee();
}

