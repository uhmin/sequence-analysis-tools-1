/* 
   translate.c
   塩基配列をアミノ酸配列に翻訳する。
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAXLETTER 256

char CODONTABLE[]="/user2/uhmin/drhamada/binsrc/universal.codon";

typedef struct{
  char *string;
  int buffer;
  int length;
}String;

typedef struct{
  FILE *fin;
  FILE *fout;
  char table[100];
  char *tablefile;
  char terminal;
  int frame;
  long from;
  long to;
  int verbose;
  int correct;
}OPTIONS;

int readfasta(FILE *fin, String *name, String *inseq);
int translate(OPTIONS options, String inseq, String *outseq);
int baseTOnum(char base);
char complement(char a);
void revcmp(char *inseq, char *revcmpseq);
int maketable(OPTIONS *options);
OPTIONS interface(int argc, char *argv[]);
int showtable(char *table);
void help(void);
String constructString(void);
int initializeString(String *data);
int deconstructString(String *data);
String fgetl(FILE *fin);

int main(int argc, char *argv[]){
  OPTIONS options;
  char inchar='\0';
  int noncode;
  String name, inseq, outseq;
  
  /* initialise */
  name=constructString();
  inseq=constructString();
  outseq=constructString();

  options=interface(argc,argv);
  maketable(&options);
  if (options.verbose>=3) showtable(options.table);

  /* rewinding */
  while(inchar!='>' && inchar!=EOF) inchar=fgetc(options.fin);
  while(feof(options.fin)==0){
    readfasta(options.fin, &name, &inseq);
    if (options.verbose>=2) fprintf(stderr,">%s%s\n",name.string, inseq.string);
    noncode=translate(options, inseq, &outseq);
    if (noncode==1 && options.verbose>=1) {
      fprintf(stderr,"\n ***** CAUTION *****\n");
      fprintf(stderr,"Found non-code codon ");
      fprintf(stderr, "(possible explanation: undecided nucleotide like N).\n");
    }
    fprintf(options.fout,">%s",name.string);
    fprintf(options.fout,"%s\n",outseq.string);
  }

  fclose (options.fin);
  fclose (options.fout);
  deconstructString(&inseq);
  deconstructString(&outseq);
  return 0;
}

String constructString(void){
  String instance;
  instance.buffer=MAXLETTER;
  instance.string=(char*)malloc(instance.buffer);
  if(instance.string==NULL){
    fprintf(stderr, "Memory err occurred on line %d.\n", __LINE__);
    exit(1);
  }
  initializeString(&instance);
  return instance;
}

int initializeString(String *data){
  if(data->string==NULL){
    data->buffer=MAXLETTER;
    data->string=(char*)malloc(data->buffer);
    if(data->string==NULL){
      fprintf(stderr, "Memory err occurred on line %d.\n", __LINE__);
      exit(1);
    }
  }
  data->string[0]='\0';
  data->length=0;
  return 0;
}

int deconstructString(String *data){
  if(data->string!=NULL){
    free(data->string);
    data->string=NULL;
  }
  data->buffer=0;
  data->length=0;
  return 0;
}

int readfasta(FILE *fin, String *name, String *inseq){
  char c;

  /* get seq name */
  initializeString(inseq);
  initializeString(name);
  *name=fgetl(fin);

  /* get sequence */
  for(c='\0'; c!='>' && c!=EOF; c=fgetc(fin)){
    if (isalpha(c)){
      if(inseq->length + 2 > inseq->buffer){
	inseq->buffer = inseq->length + MAXLETTER;
	inseq->string=(char*)realloc(inseq->string, inseq->buffer);
	if(inseq->string==NULL){
	  fprintf(stderr, "Memory err occuerred in line %d.\n", __LINE__);
	  exit(0);
	}
      }
      inseq->string[inseq->length++]=c;
    }
  }
  inseq->string[inseq->length]='\0';
  return 0;
}

int translate(OPTIONS options, String inseq, String *outseq){
  long position;
  char *tmp, base;
  int m, code, noncode=0;

  /* initialise */
  inseq.length=strlen(inseq.string);
  if (options.to==0) options.to=inseq.length;
  initializeString(outseq);

  /* frame set up */
  if (options.frame<0) {
    tmp=(char *)malloc(inseq.buffer);
    revcmp(inseq.string, tmp);
    strcpy(inseq.string, tmp);
    options.frame*=-1;
  }
  for(position=options.from-1; position < options.to; position+=3){
    code=0;
    for(m=0; m<3; m++){
      base=inseq.string[position + m + options.frame - 1];
      if(base=='\0'){
	code=-1;
	break;
      }
      code*=4;
      code+=baseTOnum(base);
    }
    if(outseq->length + 2 > outseq->buffer){
      outseq->buffer = outseq->length + MAXLETTER;
      outseq->string = (char*)realloc(outseq->string, outseq->buffer);
    }
    /* printf("%2d ",code+1); */
    if(code >= 0 && code < 64){
      outseq->string[outseq->length]=options.table[code];
    }else if(code==-1){
      break;
    }else{
      noncode=1;
      if (options.correct==1){
	if(outseq->length + 5 > outseq->buffer){
	  outseq->buffer = outseq->length + MAXLETTER;
	  outseq->string=(char *)realloc(outseq->string, outseq->buffer);
	}
        outseq->string[outseq->length++]='(';
        outseq->string[outseq->length++]=inseq.string[position + options.frame - 1];
        if (inseq.string[position + 1 + options.frame - 1] != '\0')
          outseq->string[outseq->length++] = inseq.string[position + 1 + options.frame - 1];
        if (inseq.string[position + 2 + options.frame - 1] != '\0')
          outseq->string[outseq->length++] = inseq.string[position + 2 + options.frame - 1];
        outseq->string[outseq->length]=')';
      }
      else outseq->string[outseq->length]='X';
    }
    if(outseq->string[outseq->length]=='\0') break;
    outseq->length++;
  }
  outseq->string[outseq->length]='\0';

  return noncode;
}

int baseTOnum(char base){
  base=toupper(base);
  
  switch (base){
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
  case 'U':
    return 3;
  case 0:
    return 'A';
  case 1:
    return 'C';
  case 2:
    return 'G';
  case 3:
    return 'T';
  default:
    return 100;
  }
}



char complement(char a){
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
  
  return a;
}

void revcmp(char *inseq, char *revcmpseq){
  long length,p;

  length=strlen(inseq);
  for (p=0;p<=length-1;p++) revcmpseq[p]=complement(inseq[length-p-1]);
  revcmpseq[length]='\0';
}

int maketable(OPTIONS *options){
  int code,m;
  char inchar='\0';
  FILE *fin;

  /* file open */
  if ((fin=fopen(options->tablefile,"r"))==NULL){
    fprintf(stderr,"Cannot open codonfile.[%s]\n", options->tablefile);
    exit(0);
  }
  
  while(inchar!=EOF){
    /* codon */
    code=0;
    for(m=0;m<3;m++){
      inchar=fgetc(fin);
      code*=4;
      code+=baseTOnum(inchar);
    }
    
    while(inchar!='\t' && inchar!=EOF) inchar=fgetc(fin);
    
    /* amino */
    inchar=toupper(fgetc(fin));
    if(inchar==' '){
      inchar=options->terminal;
    }
    if (code<64) {
      options->table[code]=inchar;
    }
    
    while(inchar!='\n' && inchar!=EOF) inchar=fgetc(fin);
  }
  fclose (fin);

  return 0;
}

OPTIONS interface(int argc, char *argv[]){
  int count,infile=-1,outfile=-1,tablefile=-1,region=-1;
  int n,p;
  char tmp[MAXLETTER];
  OPTIONS result;

  /* default */
  result.from=1;
  result.to=0;
  result.frame=1;
  result.tablefile=CODONTABLE;
  result.fin=stdin;
  result.fout=stdout;
  result.verbose=1;
  result.correct=0;
  result.terminal=' ';

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
    case 'f':
      result.frame=atoi(argv[count+1]);
      break;
    case 'r':
      region=count+1;
      break;
    case 't':
      tablefile=count+1;
      break;
    case 'v':
      result.verbose=atoi(argv[count+1]);
      break;
    case 'c':
      result.correct=atoi(argv[count+1]);
      if (result.correct>1 || result.correct<0) help();
      break;
    case 'T':
      result.terminal=argv[count+1][0];
      if(result.terminal=='\0'){
	result.terminal=' ';
      }
      break;
    default:
      help();
      break;
    }

  
  if (infile!=-1) {
    if ((result.fin=fopen(argv[infile],"r"))==NULL){
      fprintf(stderr,"Cannot open infile.[%s]\n",argv[infile]);
      exit(0);
    }
  }
  
  if (outfile!=-1){
    if ((result.fout=fopen(argv[outfile],"w"))==NULL){
      fprintf(stderr,"Cannot open outfile.[%s]\n",argv[outfile]);
      exit(0);
    }
  }
  if (tablefile!=-1) {
    result.tablefile=argv[tablefile];
    /* strcpy(result.table,argv[tablefile]);*/
  }

  
  /* region */
  if (region!=-1){
    n=0;
    for(p=0;argv[region][p]!='-';p++){
      tmp[n++]=argv[region][p];
      if (argv[region][p]=='\0') help();
    }
    tmp[n]='\0';
    result.from=atol(tmp);

    n=0;p++;
    for(;argv[region][p]!='\0';p++){
      tmp[n++]=argv[region][p];
    }
    tmp[n]='\0';
    result.to=atol(tmp);
  }
  
  return result;
}

int showtable(char *table){
  int one,two,three;

  fprintf(stderr,"  ---  codon table used  ---\n");
  for (one=0;one<4;one++)
    for (two=0;two<4;two++){
      for (three=0;three<4;three++)
        fprintf(stderr,"%c%c%c(%d)--%c  "
		,baseTOnum(one),baseTOnum(two),baseTOnum(three), one*16+two*4+three
               ,table[one*16+two*4+three]);
      fprintf(stderr,"\n");
    }
  printf("\n");
  return 0;
}

String fgetl(FILE *fin){
  char c;
  String result;

  /* initialize */
  result.buffer=MAXLETTER;
  result.length=0;
  result.string=(char*)malloc(result.buffer);

  /* read one line */
  c=fgetc(fin);
  for(result.length=0; c!='\n' && c!='\0' && c!=EOF; result.length++){
    result.string[result.length]=c;

    if(result.length + 4 > result.buffer){
      result.buffer = (result.length + MAXLETTER)*sizeof(char);
      result.string=(char *)realloc(result.string, result.buffer);
    }
    c=fgetc(fin);
  }

  if(c=='\n'){
    result.string[result.length] = c;
    result.length++;
  }
  result.string[result.length] = '\0';

  return result;
}

void help(void){
  fprintf(stderr,"%s written by uhmin compiled on %s %s.\n", __FILE__, __DATE__, __TIME__);
  fprintf(stderr,"function :   Translate nucleotide to amino acid.\n");
  fprintf(stderr,"  usage\n");
  fprintf(stderr,"          -i : infile, stdin if default.\n");
  fprintf(stderr,"          -o : outfile, stdout if default.\n");
  fprintf(stderr,"          -t : codontable-file, \"%s\" if default.\n",CODONTABLE);
  fprintf(stderr,"          -f : frame (-3 ~ 3), 1 if default\n");
  fprintf(stderr,"          -r : region to translate [from-to]. All [1-0] if default.\n");
  fprintf(stderr,"          -T : You can specify terminal letter (one letter such as '*').\n");
  fprintf(stderr,"               Whitespace ' ' if default.\n");
  fprintf(stderr,"          -v : information level, 1 if default.\n");
  fprintf(stderr,"               0: no informatioin.\n");
  fprintf(stderr,"               1: display only warning.\n");
  fprintf(stderr,"               2: show input nucleotide sequence.\n");
  fprintf(stderr,"               3: 2 + show codon table.\n");
  fprintf(stderr,"          -c : err correction, 0 if default.\n");
  fprintf(stderr,"               0: erroneous codon will be shown \'X\'.\n");
  fprintf(stderr,"               1: erroneoud codon will be shown like \'(NNN)\'.\n");
  fprintf(stderr,"\n");
  exit(0);
}
