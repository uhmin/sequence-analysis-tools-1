/* docono3.c find input sequence from data sequence */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#define MAXLETTER 255

void help(void);
char *mainrutin(FILE *fdata, char *indata, FILE *fout, int matchlimit, int verbose);
int compair(char *indata, char *dataword, long *position, char *prototype);
int nuccmp(char gennuc, char qnuc);
char *revcmp(char *inseq);
char revchar(char a);
char *read_indata(int infile, char *argv[]);
char *formatData(char *indata);

int main(int argc, char *argv[]){
  FILE *fout, *fdata;
  int count, datafile=0, infile=0, outfile=0, match=100,query=0, verbose=0;
  char *indata,*revdata;
  long i=0;
  char c,rev='F';
  char *result;
  
  if ((argc % 2)==0) help();
  for(count=1;count<argc-1;count+=2)
    switch(argv[count][1]){
    case 'd':
      datafile=count+1;
      break;
    case 'i':
      infile=count+1;
      break;
    case 'o':
      outfile=count+1;
      break;
    case 'm':
      match=atoi(argv[count+1]);
      break;
    case 'q':
      query=count+1;
      break;
    case 'r':
      rev=argv[count+1][0];
      break;
    case 'v':
      verbose=atoi(argv[count+1]);
      break;
    default:
      help();
    }

  if (datafile!=0){
    if ((fdata=fopen(argv[datafile],"r"))==NULL){
      fprintf(stderr,"Cannot open datafile.(%s)\n",argv[datafile]);
      exit(0);
    }
  }else {
    fprintf(stderr,"Datafile not selected.\n");
    help();
    exit(0);
  }

  if (outfile!=0){
    if ((fout=fopen(argv[outfile],"r"))==NULL){
      fprintf(stderr,"Cannot open outfile.(%s)\n",argv[outfile]);
      exit(0);
    }
  }else fout=stdout;

  if(query==0){
    indata=read_indata(infile, argv);
  }else{
    indata=(char *)malloc(sizeof(char *)*(strlen(argv[query])+2));
    for(i=0;argv[query][i]!='\0';i++)
      indata[i]=toupper(argv[query][i]);
    indata[i]='\0';
  }
  indata=formatData(indata);

  if(rev=='T') {
    /*    revdata=(char *)(sizeof(char *)*(strlen(indata)+2));*/
    revdata=revcmp(indata);
    indata=(char *)realloc(indata,sizeof(char *)*(strlen(indata)+strlen(revdata)+3));
    strcat(indata,",");
    strcat(indata,revdata);
  }
  /*  printf("%s\n",indata);*/

  while((c=fgetc(fdata))!='>') if (c==EOF) return 0;
  /*  printf("verbose=%d\n",verbose);*/
  while((result=mainrutin(fdata,indata,fout,match, verbose))!=NULL) {
    if(result[0]!='\0')
      fprintf(fout,"%s\n",result);
    free(result);
  }

  return 0;
}

char *mainrutin(FILE *fdata, char *indata, FILE *fout, int matchlimit, int verbose){
  int i,j, all_N;
  char *dataword,*prototype;
  int maxlength=0;
  long prototype_length;
  long position,preposition,qposition;
  char c;
  char *result, subresult[MAXLETTER];
  long result_position, result_length=MAXLETTER;
  int exist_result, exist_N=0;

  if(verbose==3) exist_result=1; else exist_result=0;
  /**************/
  /* initialise */
  /**************/
  /* max length of query */
  j=0;
  for(i=0;i<strlen(indata);i++){
    if (indata[i]==',' || indata[i]=='\n') {
      if (j>maxlength) maxlength=j;
      j=-1;
    }
    j++;
  }
  if (j>maxlength) maxlength=j;
  
  dataword=(char *)calloc(maxlength+5,sizeof(char *));
  prototype=(char *)calloc(maxlength+5,sizeof(char *));


  /* search sequence name */
  result=(char *)malloc(sizeof(char)*result_length);
  result_position=0;
  result[result_position++]='>';
  result[result_position]='\0';
  /*  fprintf(fout,">");  */
  while(1){
    c=fgetc(fdata);
    if (c==EOF) return NULL;
    result[result_position++]=c;
    result[result_position]='\0';
    if(result_position>=result_length){
      result_length+=MAXLETTER;
      result=realloc(result, sizeof(char *)*result_length);
    }
    if (c=='\n') break;
  }
  
  /* read first data */
  for(i=0;i<maxlength;){
    c=toupper(fgetc(fdata));
    if(c=='>') break;
    if (c!='\n') {
      dataword[i]=c;
      i++;
    }
  }
  dataword[i]='\0';
  position=1;
  preposition=0;

  /**************/
  /* main rutin */
  /**************/

  /* search and store hits */
  while(1){
    qposition=0;
    while(qposition>=0) {
      i=compair(indata,dataword,&qposition,prototype);
      if (i>=matchlimit) {
        if (preposition==0) preposition=position;
        sprintf(subresult,"%ld\t%ld\t%s\t%d%% (",position,position-preposition,prototype,i);
	result_position=strlen(subresult);
	if(verbose<=1) all_N=1; else all_N=0;
	prototype_length=strlen(prototype);
        for (j=0, exist_N=0; j<prototype_length; j++) {
	  if(dataword[j]=='N') exist_N=1; else all_N=0;
	  if(isalpha(dataword[j])){
	    subresult[result_position++]=dataword[j];
	  }
	}
	subresult[result_position]='\0';
	strcat(subresult,")\n");
        preposition=position;
	if((verbose==0 && exist_N==0) || (verbose==1 && all_N==0) || verbose>=2){
	  if(strlen(result)+strlen(subresult)>=result_length){
	    result_length+=MAXLETTER;
	    result=realloc(result, sizeof(char *)*result_length);
	  }
	  strcat(result, subresult);
	  exist_result=1;
	}
      }
    }
    for (i=0;i<maxlength;i++){
      if(dataword[i+1]=='\0') break;
      dataword[i]=dataword[i+1];
    }
    /* printf("maxlength=%d, position=%d\n", maxlength, i);*/
    while((c=fgetc(fdata))=='\n');
    dataword[i]=toupper(c);
    position++;
    if (dataword[0]==EOF || dataword[i]=='>') break;
  }

  if (dataword[0]==EOF) i=0;else i=1;
  free (dataword);
  free (prototype);
  if(!exist_result){
    result[0]='\0';
  }
  return result;
}

int compair(char *indata, char *dataword, long *position, char *prototype){
  int j,score=0;
  long length;

  length=strlen(indata);
  j=0;
  for (;*position<length;(*position)++){
    score+=nuccmp(dataword[j], indata[*position]);
    prototype[j]=indata[*position];
    if (indata[*position]==',' || indata[*position]=='\n' || indata[*position]=='\0') {
      *position=*position+1;
      break;
    }
    j++;
  }

  prototype[j]='\0';
  
  if(j==0) return 0;
  if (*position>=length) *position=-1;
  return (score*100/j);
}


int nuccmp(char gennuc, char qnuc){
  if(!isalpha(gennuc) || !isalpha(qnuc)) return 0;
  if (gennuc==qnuc) return 1;
  if (gennuc=='N' || qnuc=='N') return 1;
  
  if (qnuc=='Y' && (gennuc=='C' || gennuc=='T')) return 1;
  if (qnuc=='R' && (gennuc=='A' || gennuc=='G')) return 1;
  if (qnuc=='M' && (gennuc=='A' || gennuc=='C')) return 1;
  if (qnuc=='K' && (gennuc=='G' || gennuc=='T')) return 1;
  if (qnuc=='S' && (gennuc=='G' || gennuc=='C')) return 1;
  if (qnuc=='W' && (gennuc=='A' || gennuc=='T')) return 1;
  
  if (qnuc=='H' && (gennuc=='A' || gennuc=='C' || gennuc=='T')) return 1;
  if (qnuc=='B' && (gennuc=='G' || gennuc=='C' || gennuc=='T')) return 1;
  if (qnuc=='V' && (gennuc=='A' || gennuc=='C' || gennuc=='G')) return 1;
  if (qnuc=='D' && (gennuc=='A' || gennuc=='G' || gennuc=='T')) return 1;

  if (gennuc=='Y' && (qnuc=='C' || qnuc=='T')) return 1;
  if (gennuc=='R' && (qnuc=='A' || qnuc=='G')) return 1;
  if (gennuc=='M' && (qnuc=='A' || qnuc=='C')) return 1;
  if (gennuc=='K' && (qnuc=='G' || qnuc=='T')) return 1;
  if (gennuc=='S' && (qnuc=='G' || qnuc=='C')) return 1;
  if (gennuc=='W' && (qnuc=='A' || qnuc=='T')) return 1;

  if (gennuc=='H' && (qnuc=='A' || qnuc=='C' || qnuc=='T')) return 1;
  if (gennuc=='B' && (qnuc=='G' || qnuc=='C' || qnuc=='T')) return 1;
  if (gennuc=='V' && (qnuc=='A' || qnuc=='C' || qnuc=='G')) return 1;
  if (gennuc=='D' && (qnuc=='A' || qnuc=='G' || qnuc=='T')) return 1;  
  return 0;
}


char revchar(char a){
  if (a=='A') return ('T');
  if (a=='T') return ('A');
  if (a=='G') return ('C');
  if (a=='C') return ('G');

  if (a=='B') return ('V');
  if (a=='V') return ('B');
  if (a=='D') return ('H');
  if (a=='H') return ('D');
  if (a=='K') return ('M');
  if (a=='M') return ('K');
  if (a=='S') return ('S');
  if (a=='W') return ('W');
  if (a=='R') return ('Y');
  if (a=='Y') return ('R');

  return a;
}

char *revcmp(char *inseq){
  long length,p;
  char *outseq;
  
  /*  printf("'%s'\n",inseq);*/
  length=strlen(inseq);
  for(;length>0;length--){
    if(isalpha(inseq[length-1]) || inseq[length-1]==',') break;
  }
  outseq=(char *)malloc(sizeof(char *)*length);
  for (p=0;p<=length-1;p++) {
    outseq[p]=revchar(inseq[length-p-1]);
  }
  outseq[p]='\0';
  /*  printf("'%s'\n",outseq);exit(0);*/
  return outseq;
}

char *read_indata(int infile, char *argv[]){
  FILE *fin;
  char *indata;
  int i;
  int length=MAXLETTER;

  /* indata read */
  if (infile!=0){
    if ((fin=fopen(argv[infile],"r"))==NULL){
      fprintf(stderr,"Cannot open infile.(%s)\n",argv[infile]);
      exit(0);
    }
  }else fin=stdin;
  
  indata=(char *)malloc(length);
  for(i=0; !feof(fin);i++){
    if(i+3>length){
      length+=MAXLETTER;
      indata=(char *)realloc(indata, length);
    }
    indata[i]=toupper(fgetc(fin));
    indata[i+1]='\0';
  }
  
  fclose(fin);
  return indata;
}

char *formatData(char *indata){
  int i, j;
  int state;
  char *outdata;
  int length;

  length=strlen(indata)+10;
  outdata=(char *)malloc(length);
  for(i=0, j=0, state=0; indata[i]!='\0'; i++){
    if(j+3>length){
      length+=MAXLETTER;
      outdata=(char *)realloc(outdata, length);
    }
    if(isalpha(indata[i])){
      state=0;
      outdata[j]=indata[i];
      outdata[j+1]='\0';
      j++;
    }else if(state==0){
      strcat(outdata, ",");
      state=1;
      j++;
    }
  }
  for(;j>0;j--){
    if(isalpha(outdata[j])){
      break;
    }else{
      outdata[j]='\0';
    }
  }
  return outdata;
}

void help(void){
  printf("docono3.c written by uhmin compiled on %s\n",__DATE__);
  printf("sentence alignment tool\n");
  printf("arguments:\n");
  printf("    -d : datafile name essential.\n");
  printf("    -i : infile name stdin if default\n");
  printf("    -q : query sequence direct input\n");
  printf("    -o : outfile name stdout if default\n");
  printf("    -m : match percent (1-100) 100 if default\n");
  printf("    -r : consider reverse T/F F if default\n");
  printf("    -v : verbose 0-3, 0 if default. \n");
  printf("                0: no output if include N or no hit,\n");
  printf("                1: output if the result is not all N.\n");
  printf("                2: output even if found all N,\n");
  printf("                3: output all result anyway.\n");
  exit(0);
}
  
