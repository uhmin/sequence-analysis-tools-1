/*
  revcomp.c
  入力配列のreverse,complementを出力する。
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#define MAXLETTER 256

char complement(int a);
int interface(int argc, char *argv[], int *verbose, FILE **fin, FILE **fout, int *reverse
              , int *complement, int *width);
void help(void);
char *read_seqname(FILE *fin, char *seqname);
char *read_sequence(FILE *fin, char *sequence);
char *make_complement_sequence_of(char *sequence, long unsigned length);
char *make_reverse_sequence_of(char *sequence, long unsigned length);
int output(FILE *fout, char *seqname, char *sequence, int verbose, int width);

int main(int argc, char *argv[]){
  FILE *fin,*fout;
  char *sequence=NULL, *seqname=NULL;
  int comp,reve,verbose, width;
  long unsigned length;

  interface(argc, argv, &verbose, &fin, &fout, &reve, &comp, &width);
  while(fgetc(fin)!='>') ;   /* rewind */

  while(!feof(fin)){
    seqname=read_seqname(fin, seqname);
    sequence=read_sequence(fin, sequence);
    
    length=strlen(sequence);
    if (comp==1) make_complement_sequence_of(sequence, length);
    if (reve==1)    make_reverse_sequence_of(sequence, length);

    output(fout,seqname,sequence, verbose, width);
  }
  
  return 0;
}

char *read_seqname(FILE *fin, char *seqname){
  long i, buffsize=MAXLETTER;
  char c;
  if(seqname!=NULL) free(seqname);
  seqname=(char *)malloc(buffsize*sizeof(char*));
  for(i=0;(c=fgetc(fin))!='\n' && c!=EOF; i++){
    if(i>=buffsize){
      buffsize+=MAXLETTER;
      seqname=(char *)realloc(seqname, buffsize*sizeof(char*));
    }
    seqname[i]=c;
  }
  seqname[i]='\0';
  return seqname;
}

char *read_sequence(FILE *fin, char *sequence){
  long i, buffsize=MAXLETTER, c;
  if(sequence!=NULL) free(sequence);
  sequence=(char *)malloc(buffsize*sizeof(char*));
  for(i=0;(c=fgetc(fin))!='>' && c!=EOF;){
    if(i>=buffsize){
      buffsize+=MAXLETTER;
      sequence=(char *)realloc(sequence, buffsize*sizeof(char*));
    }
    if(isalnum(c)) {
      sequence[i]=c;
      i++;
    }
  }
  sequence[i]='\0';
  return sequence;
}

char complement(int a){
  int result, upper;
  if(isupper(a)) upper=1;
  a=tolower(a);
  switch (a){
  case 'a':
    result='t';
    break;
  case 't':
    result='a';
    break;
  case 'g':
    result='c';
    break;
  case 'c':
    result='g';
    break;
  case 'b':
    result='v';
    break;
  case 'v':
    result='b';
    break;
  case 'd':
    result='h';
    break;
  case 'h':
    result='d';
    break;
  case 'k':
    result='m';
    break;
  case 'm':
    result='k';
    break;
  case 's':
    result='w';
    break;
  case 'w':
    result='s';
    break;
  case 'r':
    result='y';
    break;
  case 'y':
    result='r';
    break;
  default:
    result='n';
  }
  if(upper) result=toupper(result);
  return result;
}

int interface(int argc, char *argv[], int *verbose, FILE **fin, FILE **fout, int *reverse
              , int *complement, int *width){
  int count,infile=-1,outfile=-1;

  /* default */
  *verbose=1;
  *reverse=1;  /* reverse = True */
  *complement=1;
  *fin=stdin;
  *fout=stdout;
  *width=60;

  /* assign */
  if ((argc % 2)==0) help();
  
  for(count=1;count<argc-1;count+=2)
    switch(argv[count][1]){
    case 'i':
      infile=count+1;
      break;
    case 'o':
      outfile=count+1;
      break;
    case 'v':
      if (strcmp(argv[count+1],"T")==0) *verbose=1;
      else if (strcmp(argv[count+1],"F")==0) *verbose=0;
      else help();
      break;
    case 'r':
      if (strcmp(argv[count+1],"T")==0) *reverse=1;
      else if (strcmp(argv[count+1],"F")==0) *reverse=0;
      else help();
      break;
    case 'c':
      if (strcmp(argv[count+1],"T")==0) *complement=1;
      else if (strcmp(argv[count+1],"F")==0) *complement=0;
      else help();
      break;
    case 'w':
      *width=atoi(argv[count+1]);
      if(*width<=0) help();
      break;
    default:
      help();
      break;
    }

  if (infile!=-1) {
    if ((*fin=fopen(argv[infile],"r"))==NULL){
      fprintf(stderr,"Cannot open infile.[%s]\n",argv[infile]);
      exit(0);
    }
  }
  
  if (outfile!=-1){
    if ((*fout=fopen(argv[outfile],"w"))==NULL){
      fprintf(stderr,"Cannot open outfile.[%s]\n",argv[outfile]);
      exit(0);
    }
  }
  
  return 0;
}

char *make_complement_sequence_of(char *sequence, long unsigned length){
  long i;
  for(i=0;i<length;i++){
    sequence[i]=complement((int)sequence[i]);
  }  
  return sequence;
}

char *make_reverse_sequence_of(char *sequence, long unsigned length){
  long i;
  char buffer;
  length--;
  for(i=0;i<length-i;i++){
    buffer=sequence[i];
    sequence[i]=sequence[length-i];
    sequence[length-i]=buffer;
  }
  return sequence;
}

int output(FILE *fout, char *seqname, char *sequence, int verbose, int width){
  int i;
  if(verbose==1){
    fprintf(stderr,"sequence length = %ld\n",strlen(sequence));
  }
  fprintf(fout,">%s",seqname);
  for(i=0; sequence[i]!='\0';i++){
    if(i % width ==0) fprintf(fout, "\n");
    fprintf(fout,"%c",sequence[i]);
  }
  fprintf(fout, "\n");
  return 0;
}

void help(void){
  fprintf(stderr,"%s written by uhmin, compiled on %s  \n", __FILE__, __DATE__);
  fprintf(stderr,"function: Reverse and/or complement input sequence\n");
  fprintf(stderr,"usage:\n");
  fprintf(stderr,"      -i: infile  name stdin  if default\n");
  fprintf(stderr,"      -o: outfile name stdout if default\n");
  fprintf(stderr,"      -r: reverse    sequence  T/F True if default\n");
  fprintf(stderr,"      -c: complement sequence  T/F True if default\n");
  fprintf(stderr,"      -v: show        message  T/F True if default\n");
  exit(0);
}

