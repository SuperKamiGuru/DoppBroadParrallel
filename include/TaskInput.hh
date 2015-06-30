#ifndef TASKINPUT_HH
#define TASKINPUT_HH

#include "MyString.h"
#include <string.h>

//MSH_include_begin
#include "MarshaledMyString.h"
//MSH_include_end

//MSH_BEGIN
class TaskInput
{
    public:
        TaskInput()
        {
            prevTemp=0.;
            newTemp=0.;
            ascii=0;
            log=0;
            overWrite=0;
            index=0;
            inDirName="";
            outDirName="";
            fileName="";
        }

        virtual ~TaskInput()
        {

        }
        TaskInput(const char *in1, const char *in2, const char *in3, double in4, double in5,
                    bool in6, bool in7, bool in8, int in9)
        {
            inDirName = *(new MyString(in1));
            outDirName = *(new MyString(in2));
            fileName = *(new MyString(in3));
            prevTemp=in4;
            newTemp=in5;
            ascii=in6;
            log=in7;
            overWrite=in8;
            index=in9;
        }
        TaskInput(const char *in1, const char *in2, const char *in3, double in4, double in5, bool in6, bool in7, bool in8)
        : index(0)
        {
            inDirName = *(new MyString(in1));
            outDirName = *(new MyString(in2));
            fileName = *(new MyString(in3));
            prevTemp=in4;
            newTemp=in5;
            ascii=in6;
            log=in7;
            overWrite=in8;
        }

        TaskInput& operator=(const TaskInput& other)
        {

            inDirName = *(new MyString(other.inDirName));
            outDirName = *(new MyString(other.outDirName));
            fileName = *(new MyString(other.fileName));
            prevTemp=other.prevTemp;
            newTemp=other.newTemp;
            ascii=other.ascii;
            log=other.log;
            overWrite=other.overWrite;
            index=other.index;

            return *this;
        }

        MyString inDirName; /* MSH: predefined */
        MyString outDirName; /* MSH: predefined */
        MyString fileName; /* MSH: predefined */

        double prevTemp; //MSH: primitive
        double newTemp; //MSH: primitive
        bool ascii; //MSH: primitive
        bool log; //MSH: primitive
        bool overWrite; //MSH: primitive
        int index; //MSH: primitive
};
//MSH_END

#endif // TASK_INPUT_HH
