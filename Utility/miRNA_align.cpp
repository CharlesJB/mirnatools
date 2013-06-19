//author: Sébastien Boisvert
// 2011-03-02
// microRNA analysis

#include <iostream>
#include <string>
#include <assert.h>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include <stdlib.h>
using namespace std;

int isPaired(string*a,int i,int j){
        char p=a->at(i);
        char q=a->at(j);
        if(p==q){
                return 0;
        }
        if(p>q){
                char t=p;
                p=q;
                q=t;
        }
        if(p=='A'&&q=='U')
                return 1;
        if(p=='C'&&q=='G')
                return 1;
        return 0;
}

string getColor(char a){
        switch (a){
                case 'A':
                        return "green";
                case 'U':
                        return "red";
                case 'C':
                        return "blue";
                case 'G':
                        return "black";
        }
        return "NULL";
}

string addColor(string a){
        ostringstream b;
        for(int i=0;i<(int)a.length();i++){
                b<<"<span style=\"color: "<<getColor(a[i])<<";\">"<<a[i]<<"</span>";
        }
        return b.str();
}

// http://www.ibluemojo.com/school/rna_folding.html
int findMaxPairs(string*a,int*B,int*P){
        int n=a->length();
        for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                        B[i*n+j]=0;
                        P[i*n+j]=0;
                }
        }
        int max;
        int tmp;
        int pos;
        for(int i=n-1;i>=0;i--){
                for(int j=0;j<n;j++){
                        max=0;
                        if(i>=j-4){
                                P[i*n+j]=-3;
                        }else{
                                tmp=B[(i+1)*n+j-1]+isPaired(a,i,j);
                                if(tmp>max){
                                        max=tmp;
                                        P[i*n+j]=-1-isPaired(a,i,j);
                                }
                                for(int k=i;k<j;k++){
                                        tmp=B[i*n+k]+B[(k+1)*n+j];
                                        if(tmp>max){
                                                max=tmp;
                                                P[i*n+j]=k;
                                        }
                                }
                                B[i*n+j]=max;
                        }
                }
        }
        return B[0*n+n-1];
}

void traceBack(int i,int j,int*B,int*P,int n,char*trace,int*p){
        if(P[i*n+j]==-3){
                for(int d=1;d<=j-i+1;d++){
                        trace[(*p)++]='.';
                }
        }else if(P[i*n+j]==-2){
                trace[(*p)++]='(';
                traceBack(i+1,j-1,B,P,n,trace,p);
                trace[(*p)++]=')';
        }else if(P[i*n+j]==-1){
                trace[(*p)++]='.';
                traceBack(i+1,j-1,B,P,n,trace,p);
                trace[(*p)++]='.';
        }else{
                int k=P[i*n+j];
                if(i<=k&&k+1<=j){
                        traceBack(i,k,B,P,n,trace,p);
                        traceBack(k+1,j,B,P,n,trace,p);
                }
        }
}

string fold(string*a){
        int B[1000000];
        int P[1000000];
        int n=a->length();
        findMaxPairs(a,B,P);
/*
        for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                        cout<<"\t"<<B[i*n+j];
                }
                cout<<endl;
        }
*/
        int maxLength=1000;
        assert(a->length()<=maxLength);
        int p=0;
        char trace[1000];
        traceBack(0,n-1,B,P,n,trace,&p);
        trace[n]='\0';
        return trace;
}

class Hit{
        int m_query;
        double m_identity;
        int m_length;
        int m_mismatches;
        int m_gaps;
        int m_queryStart;
        int m_queryEnd;
        int m_subjectStart;
        int m_subjectEnd;
public:
        void constructor(int query,double identity,int length,int mismatches,int gaps,int qStart,int qEnd,int subjectStart,
                int subjectEnd){
                m_query=query;
                m_identity=identity;
                m_length=length;
                m_mismatches=mismatches;
                m_gaps=gaps;
                m_queryStart=qStart;
                m_queryEnd=qEnd;
                m_subjectStart=subjectStart;
                m_subjectEnd=subjectEnd;
        }
        int getQuery(){
                return m_query;
        }
        int getQueryStart(){
                return m_queryStart;
        }
        int getSubjectEnd(){
                return m_subjectEnd;
        }
        int getSubjectStart(){
                return m_subjectStart;
        }
        int getQueryEnd(){
                return m_queryEnd;
        }
};

class MicroRNA{
        string m_seq;
        vector<Hit> m_hits;
        string m_key;
public:
        void constructor(string*a,string*b){
                m_key=*a;
                m_seq=*b;
        }
        void addHit(Hit*a){
                m_hits.push_back(*a);
        }
        vector<Hit>*getHits(){
                return &m_hits;
        }
        string*getSequence(){
                return &m_seq;
        }
        string*getKey(){
                return &m_key;
        }
};

class ReadSequence{
        string m_seq;
        int m_instances;
public:
        void constructor(string*a,int b){
                m_seq=*a;
                m_instances=b;
        }
        string*getSequence(){
                return &m_seq;
        }
        int getInstances(){
                return m_instances;
        }
};



int main(int argc, char* argv[]){
/*
        string a="UAACGCCAGCGUA";
        string folded=fold(&a);
        cout<<"Seq="<<a<<endl;
        cout<<"    "<<folded<<endl;
        return 0;
*/
        cout<<"microRNA miner v.1 by Sébastien Boisvert, 2011"<<endl;
        cout<<endl;
//        string readSequences="../Analysis/Unique.fasta";
//        string hairpins="../Analysis/hairpin-homo.fasta";
//        string hits="../Analysis/Hits.txt";
        //ifstream f7("../Analysis/mature-homo.fasta");
        string readSequences(argv[1]);
        string hairpins(argv[2]);
	string file_mature(argv[3]);
        string hits(argv[4]);
	string file_pairs(argv[4]);

        cout<<"Reads="<<readSequences<<endl;
        cout<<"Hairpins="<<hairpins<<endl;
        cout<<"Hits="<<hits<<endl;

        cout<<"Loading reads from "<<readSequences<<" -- please wait."<<endl;
        ifstream f1(readSequences.c_str());
        map<int,ReadSequence> reads;
        while(!f1.eof()){
                char buffer[3000];
                f1>>buffer;

                if(buffer[0]=='>'){
                        int id=atoi(buffer+1);
                        int count;
                        f1>>count;
                        string seq;
                        f1>>seq;
                        for(int i=0;i<(int)seq.length();i++)
                                if(seq[i]=='T')
                                        seq[i]='U';
                        ReadSequence a;
                        a.constructor(&seq,count);
                        reads[id]=a;
                }
        }
        f1.close();
        cout<<reads.size()<<" Reads."<<endl;

        cout<<"Loading hairpins from "<<hairpins<<" -- please wait."<<endl;
        map<string,MicroRNA> microRNAs;
        ifstream f2(hairpins.c_str());

        /*
 >hsa-let-7a-1 MI0000060 Homo sapiens let-7a-1 stem-loop
TGGGATGAGGTAGTAGGTTGTATAGTTTTAGGGTCACACCCACCACTGGGAGATAACTATACAATCTACTGTCTTTCCTA
*/
        while(!f2.eof()){
                char buffer[3000];
                f2>>buffer;

                if(buffer[0]=='>'){
                        string id=buffer+1;
                        string key;
                        f2>>key;
                        char line[3000];
                        f2.getline(line,3000);
                        string description=line;
                        string seq;
                        f2>>seq;
                        for(int i=0;i<(int)seq.length();i++)
                                if(seq[i]=='T')
                                        seq[i]='U';
                        MicroRNA a;
                        a.constructor(&key,&seq);
                        microRNAs[id]=a;
                }
        }
        f2.close();

        cout<<microRNAs.size()<<" microRNA hairpins"<<endl;

        // load mature stuff.
        map<string,string> matureStuff;

/*
 * >hsa-miR-19b-2* MIMAT0004492 Homo sapiens miR-19b-2*
 * AGTTTTGCAGGTTTGCATTTCA
 *
 */
        ifstream f7(file_mature.c_str());
        //ifstream f7("../Analysis/mature-homo.fasta");
        while(!f7.eof()){
                char buffer[3000];
                f7>>buffer;

                if(buffer[0]=='>'){
                        string id=buffer+1;
                        string key;
                        f7>>key;
                        char line[3000];
                        f7.getline(line,3000);
                        string seq;
                        f7>>seq;
                        for(int i=0;i<(int)seq.length();i++)
                                if(seq[i]=='T')
                                        seq[i]='U';
                        matureStuff[id]=seq;
                }
        }


        f7.close();

        map<string,vector< string> > hairpinToMature;
        ifstream f4(file_pairs.c_str());
        //ifstream f4("../Analysis/pairs-homo.txt");
        while(!f4.eof()){
                string a;
                string b;
                f4>>a>>b;
                hairpinToMature[a].push_back(b);
        }
        f4.close();

        ifstream f(hits.c_str());
        // headers
        // # Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score

        cout<<"Loading hits from "<<hits<<" -- please wait."<<endl;

        // 10070   hsa-mir-548f-4  95.45   22      1       0       1       22      44      23      4e-05   36.2
        while(!f.eof()){
                int query;
                string target;
                double identity;
                int length;
                int mismatches;
                int gaps;
                int queryStart;
                int queryEnd;
                int subjectStart;
                int subjectEnd;
                string trash;
                f>>query>>target>>identity>>length>>mismatches>>gaps>>queryStart>>queryEnd>>subjectStart>>subjectEnd>>trash>>trash;
                Hit h;
                h.constructor(query,identity,length,mismatches,gaps,queryStart,queryEnd,subjectStart,subjectEnd);
                if(target==""){
                        continue;
                }
                if(microRNAs.count(target)==0){
                        cout<<"Error: target not found Target='"<<target<<"'"<<endl;
                }
                assert(microRNAs.count(target)>0);

                microRNAs[target].addHit(&h);
        }

        f.close();

       
        map<int,vector<string> > microRNAsWithMaps;
        for(map<string,MicroRNA>::iterator i=microRNAs.begin();i!=microRNAs.end();i++){
                if(i->second.getHits()->size()>0){
                        map<int,map<string,int> > mapping;
                        int total=0;
                        ostringstream fileName;
                        fileName<<"html-report/"<<i->first<<".html";
                        ofstream report(fileName.str().c_str());
//                        report<<"<a href=\"../index.html\">Go back to the menu</a>";
                        report<<"<pre>"<<endl;
                        report<<i->first<<" <a href=\"http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc="<<*(i->second.getKey())<<"\">"<<*(i->second.getKey())<<"</a>"<<endl;
                        report<<addColor(*(i->second.getSequence()))<<endl;
                        report<<fold(i->second.getSequence())<<endl;

                        // show mature microRNAs
                        vector<string>*mature=&(hairpinToMature[i->first]);
                        string*hairpinSeq=i->second.getSequence();
                        for(int j=0;j<(int)mature->size();j++){
                                int space=hairpinSeq->length()+2;
                                string name=mature->at(j);
                                assert(matureStuff.count(name)>0);
                                string seq=matureStuff[name];
                                for(int offset=0;offset<hairpinSeq->size()-seq.length();offset++){
                                        int identities=0;
                                        for(int o=0;o<(int)seq.length();o++){
                                                if(seq[o]==hairpinSeq->at(o+offset)){
                                                        identities++;
                                                }
                                        }
                                        if(seq.length()-identities<3){
                                                space-=offset;
                                                while(offset--)
                                                        report<<" ";
                                                break;
                                        }
                                }
                                report<<addColor(seq);
                                space-=seq.length();
                                while(space--)
                                        report<<" ";
                                report<<name<<endl;
                        }
                        report<<endl;

                        for(int j=0;j<(int)i->second.getHits()->size();j++){
                                Hit*hit=&(i->second.getHits()->at(j));
                                int id=hit->getQuery();
                                ReadSequence*read=&(reads[id]);
                                string readSequence=*(read->getSequence());
                                int count=read->getInstances();
                                int refStart=hit->getSubjectStart();
                                int refEnd=hit->getSubjectEnd();
                                int readOffset=hit->getQueryStart();
                                int max=hit->getQueryEnd();
                                //report<<"Raw="<<readSequence<<" QFirst="<<readOffset<<" QLast="<<max<<" SFirst="<<refStart<<" SLast="<<refEnd<<endl;
                                if(refStart>refEnd){
                                        int t=refStart;
                                        refStart=refEnd-1;
                                        refEnd=t;
                                        t=max;
                                        max=readSequence.length()-readOffset+1;
                                        readOffset=readSequence.length()-t;
                                        for(int i=0;i<(int)readSequence.length();i++){
                                                char a=readSequence[i];
                                                switch (a){
                                                        case 'A':
                                                                a='U';
                                                                break;
                                                        case 'U':
                                                                a='A';
                                                                break;
                                                        case 'C':
                                                                a='G';
                                                                break;
                                                        case 'G':
                                                                a='C';
                                                                break;
                                                }
                                                readSequence[i]=a;
                                        }
                                        for(int i=0;i<(int)readSequence.length()/2;i++){
                                                char a=readSequence[i];
                                                readSequence[i]=readSequence[readSequence.length()-1-i];
                                                readSequence[readSequence.length()-1-i]=a;
                                        }
                                }
                                if(readOffset>max){
                                        max=readSequence.length()-max;
                                        readOffset=readSequence.length()-readOffset;
                                        for(int i=0;i<(int)readSequence.length();i++){
                                                char a=readSequence[i];
                                                switch (a){
                                                        case 'A':
                                                                a='U';
                                                                break;
                                                        case 'U':
                                                                a='A';
                                                                break;
                                                        case 'C':
                                                                a='G';
                                                                break;
                                                        case 'G':
                                                                a='C';
                                                                break;
                                                }
                                                readSequence[i]=a;
                                        }
                                        for(int i=0;i<(int)readSequence.length()/2;i++){
                                                char a=readSequence[i];
                                                readSequence[i]=readSequence[readSequence.length()-1-i];
                                                readSequence[readSequence.length()-1-i]=a;
                                        }


                                }
                                int h=readOffset;
                                string theRead=readSequence.substr(h,max-h);
                                mapping[refStart][theRead]+=count;
                                //report<<readSequence<<endl;
                        }

                        map<int,map<string,int> > orderedMap;

                        for(map<int,map<string,int> >::iterator ii=mapping.begin();ii!=mapping.end();ii++){
                                int pos=ii->first;
                                for(map<string,int>::iterator jj=ii->second.begin();jj!=ii->second.end();jj++){
                                        string read=jj->first;
                                        int count=jj->second;
                                        orderedMap[count][read]=pos;

                                }
                        }

                        for(map<int,map<string,int> >::reverse_iterator ii=orderedMap.rbegin();ii!=orderedMap.rend();ii++){
                                int count=ii->first;
                                for(map<string,int>::iterator jj=ii->second.begin();jj!=ii->second.end();jj++){
                                        string read=jj->first;
                                        int pos=jj->second;
                                        int l=pos;

                                        while(l--){
                                                report<<" ";
                                        }
                                        report<<addColor(read)<<" "<<count<<endl;
                                        total+=count;
                                }
                        }

                        report<<"</pre>"<<endl;
                        cout<<"Wrote "<<fileName.str()<<endl;
                        microRNAsWithMaps[total].push_back(i->first);
                }
        }

	ofstream f9("index.html");
//        ofstream f9("html-report/index.html");
//        f9<<"<h1>Metadata</h1>"<<endl;
//        f9<<"Collaborators: Patrick Provost & Aurélie Corduan<br />"<<endl;
//        f9<<"Bioinformatician: Sébastien Boisvert<br />"<<endl;
//        f9<<"<a href=\"http://genome.ulaval.ca/corbeillab/\">Jacques Corbeil Lab</a><br />";
//        f9<<"Project: microRNA mining in human cell populations using Illumina sequencing<br />"<<endl;
//        f9<<"First Meeting: Wed, 23 Feb 2011 17:10:47 -0500<br />"<<endl;
//        f9<<"Last update: Wed Mar  2 21:16:27 EST 2011<br />";
//
//
//        f9<<"<h1>Computational method</h1>"<<endl;
//        f9<<"<h2>Data preparation</h2>";
//        f9<<" MicroRNA hairpins for Homo sapiens were extracted from <a href=\"http://www.mirbase.org/\">miRBase</a> (v. 16). 1048 microRNA hairpins were obtained.";
//        f9<<" Reads were extracted from the qseq file (39502485 entries), which contained raw data including low-quality data, to produce a fasta file (also 39502485 entries).";
//        f9<<" The reads were then grouped when identical while keeping their number. There were 12297519 unique reads.";
//        f9<<" Since the reads contained 3' <a href=\"http://seqanswers.com/forums/showthread.php?t=198&highlight=sticky\">adapters</a> from the <a href=\"http://www.illumina.com/products/small_rna_sample_prep_kit.ilmn\">Illumina Small RNA Sample Prep Kit</a>, tools such as BWA could not be readily utilised as they require adapter-free sequences.";
//
//        f9<<"<h2>Read mapping</h2>";
//        f9<<" Reads were mapped with blast using parameters '-W 4 -e 0.001'. The small word length chosen (4) allowed more reads to map on microRNAs.";
//        f9<<" It was not necessary to trim reads to get rid of adapters because blast performs local alignments.";
//        f9<<" 667380 unique reads were aligned onto human microRNA hairpins. Out of 1048 hairpins, 797 had at least one alignment."<<endl;
//        f9<<" Any read may map on one or more microRNA hairpins owing to conserved hairpin seeds.";
//
//        f9<<"<h2>Secondary structures</h2> Secondary structures of microRNA hairpins were computed with the <a href=\"http://www.ncbi.nlm.nih.gov/pubmed/6161375\">Nussinov-Jacobson algorithm</a>. The latter was <a href=\"microRNA_main.cpp\">implemented in C++</a>."<<endl;

        f9<<"<h1>microRNA hairpins</h1>"<<endl;

//        f9<<"<table border=\"1\"><caption>microRNAs with Illumina tags</caption><tbody><tr><th>microRNA (hairpin)</th><th>mirbase link</th><th>Number of tags</th></tr>"<<endl;
        for(map<int,vector<string> >::reverse_iterator i=microRNAsWithMaps.rbegin();i!=microRNAsWithMaps.rend();i++){
                for(int j=0;j<i->second.size();j++){
                        string name=i->second[j];
                        f9<<"<tr><td><a href=\"html-report/"<<name<<".html\">"<<name<<"</a></td><td>";
                        string key=*(microRNAs[name].getKey());
                        f9<<"<a href=\"http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc="<<key<<"\">"<<key<<"</a></td><td>"<<i->first<<"</td></tr>"<<endl;
                }
        }
        f9<<"</tbody></table>"<<endl;
        f9.close();

        return EXIT_SUCCESS;
}
