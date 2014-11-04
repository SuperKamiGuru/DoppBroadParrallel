#ifndef TASKOUTPUT_HH
#define TASKOUTPUT_HH

//MSH_BEGIN
class TaskOutput
{
    public:
        TaskOutput()
        {
            duration=0.;
            success=0;
        }
        TaskOutput(double in1, bool in2)
        {
            duration=in1;
            success=in2;
        }

        TaskOutput& operator=(const TaskOutput& other)
        {
            duration=other.duration;
            success=other.success;
            return *this;
        }
        double duration;//MSH: primitive
        bool success;//MSH: primitive
};
//MSH_END
#endif // TASK_INPUT_HH

