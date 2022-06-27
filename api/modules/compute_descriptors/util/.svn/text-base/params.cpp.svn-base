#include "util.h"

#include <cstring>

void Params::save(const char *filename){
  cout <<"Saving parameters in "<< filename<< " ..."<< flush;
  ofstream output(filename);
  if(!output)cout << "error opening " <<filename << endl;
  for(uint i=0;i<names.size();i++){
    output << names[i] << endl;
    char *type = strchr(names[i],'.');
    if(!strcmp(type,".int") || !strcmp(type,".float"))
      output << fvalues[i] << endl;    
    else if(!strcmp(type,".char"))
      output << cvalues[i] << endl; 
  }
  output.close();
  cout << "done"<< endl;
}

void Params::load(const char *filename){
  cout <<"Loading parameters from "<< filename<< " ..."<< flush;
  ifstream input(filename);
  if(!input){
    cout << "no file " <<filename << endl;
    return;
  }
  char *bigbuf;
  float value; 
  char *bigbuf2 = new char [512];

  while(!input.eof()){    
    bigbuf = new char [512];
    input.getline(bigbuf,512);    
    char *type = strchr(bigbuf,'.');    
    if(type!=NULL){
      names.push_back(bigbuf);
      //cout << "start "<< bigbuf<< endl;
      if(!strcmp(type,".int") || !strcmp(type,".float")){
	input >> value;
	fvalues.push_back(value);
	cvalues.push_back("null");      
	input.getline(bigbuf2,512);//read end of line
      }else if(!strcmp(type,".char")){
	char *bigbuf3 = new char [512];
	input.getline(bigbuf3,512);
	fvalues.push_back(0);
	cvalues.push_back(bigbuf3);           
      } 
      //cout << type<< " " << names[names.size()-1] << " " << fvalues[names.size()-1] <<  " " << cvalues[names.size()-1] << endl; 
    }
  }
  delete []bigbuf2;
  delete []bigbuf;
  input.close();
  cout << "done"<< endl;

}

void Params::put(const char *name, float value){
   
  int  found=0;
  for(uint i=0;i<names.size();i++){
    if(!strcmp(names[i],name)){
      found=1;
      fvalues[i]=value;
    }
  }
  if(!found){
    char *buf= new char[512];
    strcpy(buf,name);
    names.push_back(buf);
    fvalues.push_back(value);
    cvalues.push_back("null");
  }  
}
void Params::put(const char *name, const char *value){
   
  int  found=0;
  for(uint i=0;i<names.size();i++){
    if(!strcmp(names[i],name)){
      found=1;
      strcpy(cvalues[i],value);
    }
  }
  if(!found){
    char *buf= new char[512];
    strcpy(buf,name);
    names.push_back(buf);
    fvalues.push_back(0);
    char *buff= new char[512];
    strcpy(buff,value);
    cvalues.push_back(buff);
  }  
}

float Params::getValue(const char *name){
  uint i=0;
  int found=0;
  float value=0;
  while(i<names.size() && !found){
    if(!strcmp(names[i],name)){
      value=fvalues[i];
      found=1;
    }
    i++;
  }
  if(!found)cout << name << " not found in params "<< endl;
  return value;
}
char * Params::getString(const char *name){
  uint i=0;
  int found=0;
  char *value=NULL;
  while(i<names.size() && !found){
    if(!strcmp(names[i],name)){
      value=cvalues[i];
      found=1;
    }
    i++;
  }
  if(!found)cout << name << " not found in params "<< endl;
  return value;
}



void loadFileNames(const char *filein, vector<char *> &filenames){

  cout << "Loading filenames from "<< filein << "... "<<flush;
  ifstream input(filein);   
  char *bigbuf;
  do{
    bigbuf = new char [512];
    input.getline(bigbuf,512);
    filenames.push_back(bigbuf);
  }while(!input.eof());
  filenames.erase((std::vector<char*>::iterator) &filenames[(int)filenames.size()-1]);
  input.close();
  cout << filenames.size()<< " files done."<< endl;
}
