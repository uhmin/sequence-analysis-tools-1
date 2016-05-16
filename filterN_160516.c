/* 
   filterN.c
   convert BLASTALL text output to tabular format.
   This program can convert not only blastn output
   but also convert blastp, tblastx output.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAXLETTER 1024

/* int sw=0; */

int verbose=1,nohits=0;
unsigned char force=0;
unsigned long Query_number=0, Database_number=0, Alignment_number=0;
unsigned long Query_count;
char qname[MAXLETTER], dname[MAXLETTER], spname[MAXLETTER], Evalue[MAXLETTER];
char *protein_name;
int protein_name_length=1;
char global_tmp[MAXLETTER], *fileend, filename[MAXLETTER]="(stdin)";;
int identity, match_length, mismatch, gaps;
long qstart, qstop, dstart, dstop, qlength, dlength;
float score;
FILE *fin, *fout, *ferr;
char AllowSameScore;

char *interface(int argc, char *argv[]);
char *search_header(void);
int search_qname_and_other(void);
char *search_query_name(void);
long search_query_length(void);
int search_dname_and_other(void);
char *get_data_name(void);
int search_identity_and_gaps(void);
long get_data_length_and_sp_name(void);
char *getdataname(int show);
int search_sp_name(void);
char *search_score_e_value(void);
int read_one_alignment(void);
char *tinypiece(void);
int get_numbers(char *indata, long *left, long *right);
int output(void);
char *treat_protein_name(char *protein_name_local, char *global_tmp_local);
char *fgets_wrap(int call, char *s, int n, FILE *stream);
char *chomp(char *text);
void help(void);

int main(int argc, char *argv[]){
  long Database_count, Alignment_count;
  int alignment_status;
  float preAlignmentScore=0, DatabaseScore=0, preDatabaseScore=0;
  unsigned int alignmentShow=1;

  interface(argc,argv);
  protein_name=(char *)malloc((size_t)protein_name_length);

  fileend=fgets_wrap(0, global_tmp,MAXLETTER-1,fin);

  Query_count=0;
  do{
  blast_start:
    Alignment_count=0;
    Database_count=0;
    preAlignmentScore=0;
    preDatabaseScore=0;
    DatabaseScore=0;
    if(fileend==NULL) goto read_end;
    /* search query name */
    protein_name[0]='\0';
    search_qname_and_other();
    Database_count=0;
    do{
      if(fileend==NULL) goto read_end;
      /* search database name */
      if(search_dname_and_other()==1) goto blast_start;
      if(fileend==NULL) goto read_end;
      if(dname[0]=='\0' || qname[0]=='\0') {
  goto blast_start;
      }
      preDatabaseScore=DatabaseScore;
      DatabaseScore=0;
      Alignment_count=0;
      alignmentShow=1;

      do{
	if(fileend==NULL) goto read_end;
	/* read alignments */
	search_sp_name();
	search_score_e_value();
	search_identity_and_gaps();
	alignment_status=read_one_alignment();
	if(DatabaseScore<1){
	  DatabaseScore=score;
	}
	/* printf("\n %ld:  %f -> %f\n",
	   Database_count, preDatabaseScore, DatabaseScore); */
	if(Database_count==Database_number && Database_number>0){
	  if(preDatabaseScore!=DatabaseScore || AllowSameScore=='F'){
	    /* printf("break!!\n"); */
	    Database_count++;
	    break;
	  }else{
	    /* printf("###### SAME!!\n"); */
	    Database_count--;
	  }
	}
	/*
	if(alignmentShow==1){
	  printf("\n %ld:  %f -> %f\n",
		 Alignment_count, preAlignmentScore, score);
	}
	*/
	if(Alignment_count>=Alignment_number && 
	   (score!=preAlignmentScore || AllowSameScore=='F')){
	  alignmentShow=0;
	}
	if(alignmentShow==1 || Alignment_number==0){
	  if(output()!=0) Alignment_count++;
	}

	preAlignmentScore=score;
      }while(strncmp(global_tmp," Score =", 8)==0 && 
             global_tmp[0]!='>' && fileend!=NULL &&
	     alignment_status!=1);
      if(Alignment_count!=0) {
	fprintf(fout,"\n");
	Database_count++;
      }
      protein_name[0]='\0';
    }while(global_tmp[0]=='>' &&
           fileend!=NULL &&
           (Database_number==0 || Database_count<=Database_number)&&
	   alignment_status!=1);

    if(Database_count!=0) {
      fprintf(fout,"\n");
    }
    search_header();
  read_end:
    ;
  }while(fileend!=NULL && Query_count!=Query_number);


  if (verbose>=1){
    fprintf(ferr,"filename %s output %ld sequences",filename,(Query_count-nohits));
    if (nohits!=0) fprintf(ferr,", no hits %d, total %ld",nohits, Query_count);
    fprintf(ferr,".\n");
  }
  return 0;
}


char *interface(int argc, char *argv[]){
  int count,infile=0,outfile=0,errfile=0;
  AllowSameScore='F';
  force=0;

  if ((argc % 2)==0) help();
  for(count=1;count<argc-1;count+=2){
    if(strcmp(argv[count],"-i")==0){
      infile=count+1;
    }else if(strcmp(argv[count],"-o")==0){
      outfile=count+1;
    }else if(strcmp(argv[count],"-e")==0){
      errfile=count+1;
    }else if(strcmp(argv[count],"-v")==0){
      verbose=atoi(argv[count+1]);
    }else if(strcmp(argv[count],"-Q")==0){
      Query_number=atol(argv[count+1]);
    }else if(strcmp(argv[count],"-D")==0){
      Database_number=atol(argv[count+1]);
    }else if(strcmp(argv[count],"-A")==0){
      Alignment_number=atol(argv[count+1]);
    }else if(strcmp(argv[count],"-f")==0){
      if(strcmp(argv[count+1],"T")==0){
        force=1;
      }
    }else if(strcmp(argv[count],"-S")==0){
      if(strcmp(argv[count+1],"T")==0){
        AllowSameScore='T';
      }else if(strcmp(argv[count+1],"F")==0){
        AllowSameScore='F';
      }else{
	fprintf(stderr, "filterN: invalid value on option -S. Please input T or F\n");
        help();
      }
    }else{
      fprintf(stderr,"filterN: invalid option  (%s)\n", argv[count]); 
      help();
    }
  }
  if (errfile!=0){
    if ((ferr=fopen(argv[errfile],"w"))==NULL){
      fprintf(stderr,"Cannot open errfile.(%s)\n",argv[errfile]);
      exit(0);
    }
  }else ferr=stderr;
  
  if (infile!=0){
    strcpy(filename,argv[infile]);
    if ((fin=fopen(argv[infile],"r"))==NULL){
      fprintf(ferr,"Cannot open infile.(%s)\n",argv[infile]);
      exit(0);
    }
  }else fin=stdin;

  if (outfile!=0){
    if ((fout=fopen(argv[outfile],"w"))==NULL){
      fprintf(ferr,"Cannot open outfile.(%s)\n",argv[outfile]);
      exit(0);
    }
  }else fout=stdout;

  return filename;
}


int search_qname_and_other(void){
    qname[0]='\0';
    search_query_name();
    search_query_length();
    return 0;
}

char *search_header(void){
  while(fileend!=NULL){
    /*  printf("%s",global_tmp);*/
    if (strncmp(global_tmp,"Query= ", 7)==0) {
      break;
    }
    fileend=fgets_wrap(0, global_tmp,MAXLETTER-1,fin);
  }
  return fileend;
}

char *search_query_name(void){
  /*  char *cp1;*/
  char *name_search_qname;
  while(fileend!=NULL){
    if(strncmp(global_tmp,"Query= ", 7)==0) break;
    /*    if (cp1!=NULL) break;*/
    fileend=fgets_wrap(0, global_tmp,MAXLETTER-1,fin);
  }
  if(fileend==NULL) return NULL;
  Query_count++;
  /* get query name */
  /* printf("qname length = %d\n",strlen(global_tmp)); */
  name_search_qname=strtok((global_tmp+7)," \t\n");
  if(name_search_qname!=NULL){
    strcpy(qname,name_search_qname);
  }else{
    strcpy(qname, " ");
  }
  /* printf("queryname=(%s)\n",qname); */
  return qname;
}


long search_query_length(void){
  char *cp2;
  int n,m;

  while(1){
    if(strstr(global_tmp,"letters)")!=NULL){
      cp2=strtok((strstr(global_tmp,"(")+1)," ");
      break;
    }else if(strncmp(global_tmp, "Length=", 7)==0){
      cp2=strstr(global_tmp,"=")+1;
      //fprintf(stderr, "found length: %s\n", cp2);
      break;
    }
    fileend=fgets_wrap(__LINE__, global_tmp,MAXLETTER-1,fin);
    if(fileend==NULL) return 0;
  }
  
  /*
   * if order of length over 10000 "," is inserted
   * the following process will delete non-number charactor
  */
  for(n=0,m=0;n<100;n++){
    if (isdigit(cp2[n])) {
      global_tmp[m++]=cp2[n];
      global_tmp[m]='\0';
    }
    if (cp2[n]=='\0') break;
  }
  /* get query length */
  qlength=atol(global_tmp);
  /*printf("query length=%ld\n",qlength);*/
  return qlength;
}


int search_dname_and_other(void){
  char *data_name;
  dname[0]='\0';
  data_name=get_data_name();
  if(data_name==NULL) return 1;
  else if(strcmp(data_name,"NO HITS")!=0 && fileend!=NULL){
    get_data_length_and_sp_name();
  }
  return 0;
}


char *get_data_name(void){
  char *cp1;

  /* search data name */
  while(global_tmp[0]!='>'){
    if(strncmp(global_tmp, "Query= ", 7)==0) {
      return NULL;
    }
    if (strstr(global_tmp,"No hits found")!=NULL || fileend==NULL) {
      nohits++;
      if (verbose>=2)
        fprintf(ferr,"%s no blast hits\n",qname);
      return "NO HITS";
    }
    fileend=fgets_wrap(__LINE__, global_tmp,MAXLETTER-1,fin);
  }
  
  /* get data name */
  cp1=strtok((global_tmp+1)," \t\n");
  strcpy(dname,cp1);
  /* printf("dataname=%s\n",dname); */

  /* for searching species name */
  cp1=strtok(NULL,"\n");
  if (cp1!='\0') strcpy(global_tmp,cp1); else global_tmp[0]='\0';

  return dname;
}


long get_data_length_and_sp_name(void){
  char *cp1;

  spname[0]='\0';
  while(fileend!=NULL){
    /*    search_sp_name(global_tmp,spname);*/
    search_sp_name();
    if(strncmp(global_tmp, "Length=", 7)==0){
      //fprintf(stderr, "FOUND: %s\n", global_tmp);
      //fprintf(stderr, "data length: %s\n", global_tmp+7);
      dlength=atol(global_tmp+7);
      //fprintf(stderr, "int data length: %ld\n", dlength);
      break;
    }else{
      cp1=strstr(global_tmp,"Length =");
      if (cp1!=NULL){
	dlength=atol(cp1+8);
	break;
      }
    }
    fileend=fgets_wrap(__LINE__, global_tmp,MAXLETTER-1,fin);
  }
  /*printf("datalength=%ld\n",dlength);*/

  return dlength;
}


int search_sp_name(void){
  char *cp1,*cp2;
  long i=0;
  /* get species name */
  if ((cp1=strstr(global_tmp,"["))!=NULL){
    protein_name=treat_protein_name(protein_name, global_tmp);
    cp2=strtok(cp1+1,"]\n");
    strcpy(spname,cp2);
  }else if (strstr(global_tmp,"]")!=NULL){
    i=0;
    while(isspace(global_tmp[i++])) {    }
    strtok(global_tmp,"]\n");
    cp2=global_tmp+i-2;
    strcat(spname,cp2);
  }else if(strstr(global_tmp,"Length = ")==NULL
	   && strstr(global_tmp," Score = ")==NULL ){
   protein_name=treat_protein_name(protein_name, global_tmp);
  }
  return 0;
}


char *search_score_e_value(void){
  char *cp1, *cp2;

  while(fileend!=NULL){
    if(strncmp(global_tmp," Score = ", 9)==0) break;
    fileend=fgets_wrap(__LINE__, global_tmp,MAXLETTER-1,fin);
  }
  /* get score */
  score=atof(global_tmp+8);
  cp1=strstr(global_tmp,"Expect");

  /*get e-value */
  cp2=strstr(cp1,"=")+1;
  while(isspace(*cp2)) cp2++;
  cp1=strtok(cp2,",\t\n");
  strcpy(Evalue,cp1);
  
  return NULL;
}


int search_identity_and_gaps(void){
  char *cp1, *cp2;
  long match;

  /* search identities, gaps */
  fileend=fgets_wrap(__LINE__, global_tmp,MAXLETTER-1,fin);
  if(fileend==NULL) return 0;

  cp1=strstr(global_tmp,"Identities =")+12;
  /* get match */
  match=atol(strstr(global_tmp,"Identities =")+12);
  cp2=strstr(cp1,"/")+1;

  /* get identity */
  match_length=atol(cp2);
  identity=100*match/match_length;
  /*printf("identities= %ld/%ld(%d%%),  ",match,match_length,identity);*/

  mismatch=match_length-match;

  /* get gaps if exists */
  cp1=strstr(global_tmp,"Gaps =");
  if (cp1!=NULL){
    gaps=atol(cp1+6);
  }else gaps=0;
  /*printf("gaps= %ld\n",gaps);*/

  return 0;
}


int read_one_alignment(void){
  long qstart_tmp, dstart_tmp;
  int result=0;

  fileend=tinypiece();
  qstart_tmp=qstart, dstart_tmp=dstart;

  while(fileend!=NULL && *fileend!='N' && *fileend!='B')
    fileend=tinypiece();
  qstart=qstart_tmp, dstart=dstart_tmp;
  if(*fileend=='B') result=1;
  return result;
}


char *tinypiece(void){
  while(1){
    /*    fprintf(stderr, "%s", global_tmp);*/
    /* if(strstr(global_tmp, "AK084006.1")) sw=1; */
    /*    cp1=strstr(global_tmp,"Query:");*/
    /*    if (cp1!=NULL) break;*/
    if(strncmp(global_tmp,"Query:", 6)==0) break;
    if (strncmp(global_tmp,"Query= ", 7)==0
	) {
      /* fprintf(stderr, "found BLAST\n"); */
      return "B";
    }
    if (strncmp(global_tmp," Score =", 8)==0 || global_tmp[0]=='>'
        || strstr(global_tmp,"Database:")!=NULL) {
      return "N";
    }
    fileend=fgets_wrap(0, global_tmp,MAXLETTER-1,fin);
    if(fileend==NULL) return NULL;
  }

  /* get qstart and qstop */
  get_numbers(global_tmp,&qstart,&qstop);
  fileend=fgets_wrap(__LINE__, global_tmp,MAXLETTER-1,fin);
  fileend=fgets_wrap(__LINE__, global_tmp,MAXLETTER-1,fin);
  if(fileend==NULL) return NULL;
  /* get dstart and dstop */
  get_numbers(global_tmp,&dstart,&dstop);
  
  /*printf("qstart=%ld, qstop=%ld, dstart=%ld, dstop=%ld\n",*qstart,*qstop,*dstart,*dstop);*/

  return "1";
}

int get_numbers(char *indata, long *left, long *right){

  while(*indata!='\0') {
    if(isdigit(*indata)!=0) break;
    else indata++;
  }
  *left=atol(indata);
  while(*indata!='\0') {
    if(isdigit(*indata)==0) break;
    else indata++;
  }
  while(*indata!='\0') {
    if(isdigit(*indata)!=0) break;
    else indata++;
  }
  *right=atol(indata);
  /*printf("left=%ld, right=%ld\n",*left,*right);*/
  return 0;
}


int output(void){
  if(qname[0]=='\0' || dname[0]=='\0') {
    fileend=NULL;
    return 0;
  }

    fprintf(fout,"%s\t%s\t",qname,dname);
    fprintf(fout,"%d\t%d\t",identity,match_length);
    fprintf(fout,"%d\t%d\t",mismatch,gaps);
    fprintf(fout,"%ld\t%ld\t",qstart,qstop);
    fprintf(fout,"%ld\t%ld\t",dstart,dstop);
    fprintf(fout,"%s\t%.1f\t",Evalue,score);
    fprintf(fout,"%ld\t%ld\t",qlength,dlength);
    fprintf(fout,"%s\t%s\t",spname, protein_name);
    fprintf(fout,"\n");
    /*    protein_name[0]='\0';*/
    return 1;
}

char *treat_protein_name(char *protein_name_local, char *global_tmp_local){
  int i=0;
  char *cp1;

  i=0;
  while(isspace(global_tmp_local[i++])) {}
  i-=2; if(i<=0) i=0;
  cp1=global_tmp_local+i;
  if(strlen(cp1)+strlen(protein_name_local)>protein_name_length){
    protein_name_length+=MAXLETTER;
    protein_name_local=(char *)realloc(protein_name_local, (size_t)protein_name_length);
  }
  strcat(protein_name_local, cp1);
  strtok(protein_name_local,"[");
  chomp(protein_name_local);
  return protein_name_local;
}

char *fgets_wrap(int call, char *s, int n, FILE *stream){
  char *result;
  static unsigned long line;

  line++;
  result=fgets(s, n, stream);
  if(force==1){
    call=0;
  }
  if(call != 0 
     && (strncmp(s, "Query= ", 7)==0)){
    fprintf(stderr, "s: %s\n", s);
    fprintf(stderr, "Data err! Format of input data may be wrong near line %ld\n", line);
    fprintf(stderr, "Called at line %d\n", call);
    exit(0);
  }
  /*  fprintf(stderr, "pass\n");*/
  return result;
}


char *chomp(char *text){
  int i;
  i=strlen(text);
  for(; i>=0; i--){
    if(text[i]=='\r' || text[i]=='\n'){
      text[i]='\0';
    }
  }
  return text;
}

void help(void){
  fprintf(stderr,"%s written by Uhmin. Compiled on %s %s\n", __FILE__, __DATE__, __TIME__);
  fprintf(stderr,"function : Convert blast textbase output into tabular format.\n");
  fprintf(stderr,"swhitches\n");
  fprintf(stderr,"       -i : infile name. stdin if default.\n");
  fprintf(stderr,"       -o : outfile name. stdout if default.\n");
  fprintf(stderr,"       -e : errmessege output. stderr if default.\n");
  fprintf(stderr,"       -v : message level (0-2). 1 if default.\n");
  fprintf(stderr,"       -Q : number of query  sequence  to show. If 0 show all. 0 if default\n");
  fprintf(stderr,"       -D : number of matched database to show. If 0 show all. 0 if default\n");
  fprintf(stderr,"       -A : number of aligments to show.        If 0 show all. 0 if default\n");
  fprintf(stderr,"       -S : allow Same score [T/F]. %c if default. If select T allow more than -D specified number if the result show the Same score.\n", AllowSameScore);
  fprintf(stderr, "\n\n");
  fprintf(stderr, "***** CAUTION ******\n");
  fprintf(stderr, "This program may find BLAST format err.\n");
  fprintf(stderr, "If format error were found, this program will stop processing.\n");
  fprintf(stderr, "If you DO WANT TO FORCE processing, try '-f T' option.\n");
  fprintf(stderr, "\n");
  exit(0);
}
