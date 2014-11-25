#ifndef MyString_h
#define MyString_h

#include <string.h>

//MSH_BEGIN
class MyString{
public:
    int length() { return len;}
    char* buffer() { return str;}
    MyString() { str = NULL; len = 0;}
    MyString(const char* s){
        int i=0;
        //find the length of the string
        while(s[i]!=0) i++;
        len = i;
        str = (char*)malloc(len + 1);
        memcpy(str,s,len);
        *(str+len)='\0';
    }
    ~MyString(){
	if(str!=NULL) free(str);
    }
    MyString& operator=(MyString source)
    {
	if (this != &source) {
	    // copy all our own fields here.
	    len = source.length();
	    str = (char*)malloc(len + 1);
	    memcpy(str,source.buffer(),len);
	    *(str+len)='\0';
	}
	return *this;
    }//end operator=

    bool equals(char* s){
	if(strcmp(str,s)==0) return true;
	return false;
    }
    string GetString()
    {
        string result="";
        for(int i=0; i<len; i++)
        {
            result.push_back(str[i]);
        }
        return result;
    }

protected:
    int len;	// MSH: primitive
    char* str;  /* MSH: manual
		  {
	             memcpy($$, (char*)$THIS->buffer(), $THIS->length());
		     *($$+$THIS->length()) = '\0';
		  }
		  {
	             MyString* tmp = new MyString($$);
		     memcpy($THIS, tmp, sizeof(MyString));
		  }
		  {
		     $SIZE = $THIS->length()+1; //one byte for the null-termination
		  } */
};
//MSH_END
#endif
