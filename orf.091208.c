/*
  orf.c
  入力配列6フレームの中から最長のORFを見つけて出力する。
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#define MAXLETTER 256
#define CODONLETTER 30
#define START 0
#define STOP 1

typedef struct{
  char *string;
  long buffer;
  long length;
}STRING;

typedef struct{
  STRING name;
  STRING sequence;
  STRING orf;
  long from;
  long to;
}SEQUENCE;

typedef struct{
  int verbose;
  FILE *fin;
  FILE *fout;
  int reverse;
  char *startCodon;
  char *stopCodon;
  int addNothing;
  int terminal;
}OPTIONS;

OPTIONS interface(int argc, char *argv[]);
SEQUENCE readOneSequence(OPTIONS options);
int searchORF(SEQUENCE *SEQ ,OPTIONS options);
int ORFoneframe(char *sequence, OPTIONS options, int shift, long *start, long *stop);
STRING *copyorf(STRING *orf, char *sequence, long start, long stop);
void revcmp(char *inseq, char *outseq);
char revchar(char a);
int nuccmp(char gennuc_p, char qnuc_p);
int codonMatch(char *test, char *codon, int sw);
void help(void);
STRING instanceString(void);
STRING fgetl(FILE *fin);
void freeSEQUENCE(SEQUENCE seq);
void freeSTRING(STRING line);

int main(int argc, char *argv[]){
  SEQUENCE SEQ;
  OPTIONS options;
  int orientation;
  int from, to;

  options=interface(argc,argv);
  if (options.verbose>=3){
    fprintf(options.fout,"# number of start codon= %d\n"
	    ,(strlen(options.startCodon)+1)/4);
  }

  /* search start opint of sequence name */
  while(fgetc(options.fin)!='>' && !feof(options.fin)) ;
  while(1){
    SEQ=readOneSequence(options);
    /* main program */
    orientation=searchORF(&SEQ ,options);
    fprintf(options.fout,">%s",SEQ.name.string);

    if(SEQ.orf.string!=NULL){
      if(orientation<0) fprintf(options.fout," (complement)");
      if(options.addNothing==0){
	from=SEQ.from+1;
	to=SEQ.to+3;
	if(from>SEQ.sequence.length){
	  from=SEQ.sequence.length;
	}
	if(to>SEQ.sequence.length){
	  to=SEQ.sequence.length;
	}
	fprintf(options.fout," from=%d to=%d length=%d"
		, from, to, to-from+1);
      }
      fprintf(options.fout,"\n");
      fprintf(options.fout,"%s\n",SEQ.orf.string);
    }else{
      fprintf(options.fout,"  -- No ORF found --\n\n");
    }
    if (options.verbose>=1) fprintf(options.fout,"#\n");
    SEQ.sequence.string[0]='\0';
    if (feof(options.fin)) break;
    freeSEQUENCE(SEQ);
  }

  fclose (options.fin);
  fclose (options.fout);
  return 0;
}

SEQUENCE readOneSequence(OPTIONS options){
  SEQUENCE SEQ;
  char inchar;

  /* input sequence name */
  SEQ.name=fgetl(options.fin);
  SEQ.name.string[--SEQ.name.length]='\0'; /* chop */
  SEQ.orf.string=NULL;
  
  if (options.verbose>=1) fprintf(options.fout,"# session %s\n", SEQ.name.string);

  /* read sequence */
  inchar=fgetc(options.fin);
  SEQ.sequence=instanceString();
  for(SEQ.sequence.length=0; inchar!=EOF;){
    /*if (('a'<=inchar && inchar<='z') || ('A'<=inchar && inchar<='Z')) {*/
    if (isalpha(inchar)) {
      if(SEQ.sequence.length + 10 > SEQ.sequence.buffer){
	SEQ.sequence.buffer = (SEQ.sequence.length + MAXLETTER)*sizeof(char);
	SEQ.sequence.string 
	  = (char*)realloc(SEQ.sequence.string, SEQ.sequence.buffer);
	if(SEQ.sequence.string == NULL){
	  fprintf(stderr, "Memory err occurred on line %d\n", __LINE__);
	  exit(1);
	}
      }
      SEQ.sequence.string[SEQ.sequence.length]   = inchar;
      SEQ.sequence.string[SEQ.sequence.length+1] = '\0';
      SEQ.sequence.length++;
    }
    inchar=fgetc(options.fin);
    if (inchar=='>') break;
  }
  return SEQ;
}

int searchORF(SEQUENCE *SEQ ,OPTIONS options){
  STRING revSequence;
  long start=0, stop=0;
  int maxframe=99, frame2;
  int from, to;

  SEQ->from=0;
  SEQ->to=0;

  /* print information */
  if(options.verbose>=1){
    fprintf(options.fout,"#         length of query sequence: %ld\n", SEQ->sequence.length);
  }
  if (options.verbose>=4){
    fprintf(options.fout,"# query sequence:\n");
    fprintf(options.fout,"# %s\n",SEQ->sequence.string);
  }
  
  /* for overse sequence */
  for (frame2=0; frame2 < 3; frame2++){
    if (options.verbose>=2){
      fprintf(options.fout,"# ***   frame=%d\n",frame2+1);
    }
    ORFoneframe(SEQ->sequence.string, options, frame2, &start, &stop);
    if ((stop-start) > (SEQ->to - SEQ->from)) {
      SEQ->from = start;
      SEQ->to   = stop;
      maxframe=frame2+1;
      copyorf(&SEQ->orf, SEQ->sequence.string, start, stop);
    }
  }

  /* for reverse sequence */
  /* make reverse_complement */
  revSequence.buffer=(SEQ->sequence.length+1)*sizeof(char);
  revSequence.string=(char *)malloc(revSequence.buffer);
  revcmp(SEQ->sequence.string, revSequence.string);
  for (frame2=0; frame2 < 3*options.reverse; frame2++){
    if (options.verbose>=2){
      fprintf(options.fout,"# ***   frame=%d\n",-1-frame2);
    }
    ORFoneframe(revSequence.string, options, frame2, &start, &stop);
    if ((stop-start)>(SEQ->to - SEQ->from)) {
      SEQ->from = start;
      SEQ->to   = stop;
      maxframe=-1-(frame2);
      copyorf(&SEQ->orf, revSequence.string, start, stop);
    }
  }
  
  if (options.verbose >= 1){
    from = SEQ->from+1;
    to   = SEQ->to+3;
    if(from > SEQ->sequence.length){
      from = SEQ->sequence.length;
    }
    if(to > SEQ->sequence.length){
      to = SEQ->sequence.length;
    }
    fprintf(options.fout,
	    "# _____  selected frame= %d, from %d to %d, length %d bp  _____\n"
            ,maxframe, from, to, to-from+1);
  }
  free(revSequence.string);
  revSequence.string=NULL;
  return maxframe;
}

int ORFoneframe(char *sequence, OPTIONS options, int shift, long *start, long *stop){
  long orf[2]={0,0};
  long unsigned p;
  int orfsw=0, variety, codonnum;
  char *startcodon, *stopcodon;
  startcodon = options.startCodon;
  stopcodon  = options.stopCodon;

  codonnum=(strlen(startcodon)+1)/4;
  for(p=shift; p<strlen(sequence); p+=3){
    /* search start codon */
    for(variety=0; variety < codonnum; variety++){
      if (orfsw==0 && codonMatch(sequence+p, startcodon+variety*4, 0)==1){
	orfsw=1;
	*start=p;
	p+=3;
      }
    }
    
    codonnum=(strlen(stopcodon)+1)/4+1;
    for(variety=0; variety < codonnum; variety++){
      /* search stop codon */
      if (orfsw==1 && codonMatch(sequence+p, stopcodon+variety*4, options.terminal)==1){
	orfsw=0;
	*stop=p;
	if (options.verbose >= 3)
	  fprintf(options.fout,"#        %ld - %ld length= %ld\n"
		  ,*start+1,*stop+3,*stop-*start+3);
	if ((*stop-*start) > (orf[STOP]-orf[START])) {
	  orf[STOP]  = *stop;
	  orf[START] = *start;
	}
      }
    }
    if (sequence[p]=='\0') break;
  }

  *start=orf[START];
  *stop =orf[STOP];

  if (options.verbose>=2 && orf[START]>=0 && orf[STOP]>=0 && orf[START]!=orf[STOP])
    fprintf(options.fout,"#       query length=%d,  positioin  %ld - %ld,   orf = %ld bp\n"
            , strlen(sequence), orf[START]+1, orf[STOP]+3, orf[STOP]-orf[START]+3);
  else if (options.verbose>=2)
    fprintf(options.fout,"#       query length=%d,               no hits found\n"
            ,strlen(sequence));

  return 0;
}

STRING *copyorf(STRING *orf, char *sequence, long start, long stop){
  long n;

  if(orf->string != NULL){
    free(orf->string);
    orf->string = NULL;
  }
  orf->buffer = (stop-start+10)*sizeof(char);;
  orf->string = (char *)malloc(orf->buffer);
  if(orf->string == NULL){
    fprintf(stderr, "Memory err occurred on line %d\n", __LINE__);
    exit(1);
  }
  orf->string[0] = '\0';
  for(n=start, orf->length=0; n<stop+3; n++) {
    orf->string[orf->length++]=sequence[n];
  }
  orf->string[orf->length]='\0';

  return orf;
}


char revchar(char a){
  switch(a){
  case 'A': return('T'); break;
  case 'T': return('A'); break;
  case 'G': return('C'); break;
  case 'C': return('G'); break;
  case 'B': return('V');  break;
  case 'V': return('B'); break;
  case 'D': return('H'); break;
  case 'H': return('D'); break;
  case 'K': return('M'); break;
  case 'M': return('K'); break;
  case 'S': return('S'); break;
  case 'W': return('W'); break;
  default:  return(a);
  }
  return a;
}

void revcmp(char *inseq, char *outseq){
  long length,p;
  char c;
  
  length = strlen(inseq) -1;
  for (p=0; p<=length; ) {
    c = revchar(inseq[length-p]);
    if (isalpha(c)) outseq[p++] = c;
  }
  outseq[p] = '\0';
}

OPTIONS interface(int argc, char *argv[]){
  int count,infile=-1,outfile=-1;
  int length;
  OPTIONS result;

  /* default */
  result.verbose=1;
  result.reverse=1;  /* reverse = True */
  result.fin=stdin;
  result.fout=stdout;
  result.addNothing=0;
  result.startCodon=NULL;
  result.stopCodon=NULL;
  result.terminal=1;

  /* assign */
  if ((argc % 2)==0) help();
  
  for(count=1; count<argc-1; count+=2)
    switch(argv[count][1]){
    case 'i':
      infile=count+1;
      break;
    case 'o':
      outfile=count+1;
      break;
    case 'v':
      result.verbose=atoi(argv[count+1]);
      break;
    case 'r':
      if (strcmp(argv[count+1],"T")==0) result.reverse=1;
      else if (strcmp(argv[count+1],"F")==0) result.reverse=0;
      else help();
      break;
    case 'S':
      length=strlen(argv[count+1]);
      result.startCodon=(char*)malloc(length+1);
      strcpy(result.startCodon, argv[count+1]);
      break;
    case 'E':
      length=strlen(argv[count+1]);
      result.stopCodon=(char*)calloc(length+10, sizeof(char));
      strcpy(result.stopCodon, argv[count+1]);
      break;
    case 'a':
      if (strcmp(argv[count+1],"T")) result.addNothing=1;
      else if (strcmp(argv[count+1],"F")) result.addNothing=0;
      break;
    case 'T':
      if (strcmp(argv[count+1],"T")==0) result.terminal=1;
      else if (strcmp(argv[count+1],"F")==0) result.terminal=0;
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

  if(result.startCodon==NULL){
    result.startCodon=(char*)malloc(CODONLETTER);
    strcpy(result.startCodon,"ATG,GTG");
  }
  if(result.stopCodon==NULL){
    result.stopCodon=(char*)calloc(CODONLETTER, sizeof(char));
    strcpy(result.stopCodon,"TGA,TAG,TAA");
  }
  return result;
}

void help(void){
  fprintf(stderr,"%s written by uhmin, compiled on %s %s\n", __FILE__, __DATE__, __TIME__);
  fprintf(stderr,"function   : Select and output longest ORF from inserted sequence.\n");
  fprintf(stderr,"             Only ACGT system is available.\n");
  fprintf(stderr,"arguments  :\n");
  fprintf(stderr,"            -i: Infile  name stdin if default.\n");
  fprintf(stderr,"            -o: Outfile name stdin if default.\n");
  fprintf(stderr,"            -r: If onsider reverse sequence. [T/F] T if default.\n");
  fprintf(stderr,"            -T: Consider sequence end as terminal codon [T/F]. T if default \n");
  fprintf(stderr,"            -S: Start codon \"ATG,GTG\" if default.\n");
  fprintf(stderr,"            -E: Stop  codon \"TGA,TAG,TAA\" if default.\n");
  fprintf(stderr,"            -v: Message level to show. 1 if default.\n");
  fprintf(stderr,"                -v 0: Output only ORF sequence.\n");
  fprintf(stderr,"                -v 1: -v 0 + the longest ORF position.\n");
  fprintf(stderr,"                -v 2: -v 1 + the longest ORF position for every frame.\n");
  fprintf(stderr,"                -v 3: -v 2 + all the ORF position found.\n");
  fprintf(stderr,"                -v 4: -v 3 + query sequence.\n");
  fprintf(stderr,"            -a: Add information to sequence name [T/F], T if default.\n");
  exit(0);
}

int codonMatch(char *test, char *codon, int sw){
  int result=0;
  if(nuccmp(test[0], codon[0])==1 && 
     nuccmp(test[1], codon[1])==1 && 
     nuccmp(test[2], codon[2])==1){
    result=1;
  }else if(sw==1 && codon[0]=='\0' && 
	   (test[3]=='\0' || test[2]=='\0' || test[1]=='\0' || test[0]=='\0')){
    result=1;
  }

  return result;
}


int nuccmp(char gennuc_p, char qnuc_p){
  char gennuc, qnuc;
  gennuc=toupper(gennuc_p);
  qnuc=toupper(qnuc_p);
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

STRING instanceString(void){
  STRING data;

  data.buffer=MAXLETTER;
  data.length=0;
  data.string=(char*)malloc(data.buffer);
  if(data.string==NULL){
    fprintf(stderr, "Data err on line %d\n", __LINE__);
    exit(1);
  }
  data.string[0]='\0';
  return data;
}

STRING fgetl(FILE *fin){
  char c;
  STRING result;

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

void freeSEQUENCE(SEQUENCE seq){
  freeSTRING(seq.name);
  freeSTRING(seq.sequence);
  freeSTRING(seq.orf);
  seq.from = 0;
  seq.to   = 0;
}
void freeSTRING(STRING line){
  if(line.string!=NULL){
    free(line.string);
    line.string=NULL;
  }
  line.buffer=0;
  line.length=0;
}
