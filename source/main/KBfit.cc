#include "task_handler.h"
#include <cstdio>
#include <ctime>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

// ******************************************************************************
// * *
// *                 Main driver program to run "KBfit" *
// * *
// *   Program takes a single argument that is the name of the input file. *
// *   Input file must contain a single XML document with root tag named *
// *   <KBFit>.  The input XML must have the form below: *
// * *
// *    <KBFit> *
// * *
// *       <Initialize> *
// *         <ProjectName>NameOfProject</ProjectName> *
// *         <LogFile>output.log</LogFile> *
// *         <EchoXML/> *
// *         <MCSamplingInfo> ... </MCSamplingInfo> *
// *       </Initialize> *
// * *
// *       <TaskSequence> *
// *         <Task><Action>...</Action> ...  </Task> *
// *         <Task><Action>...</Action> ...  </Task> *
// *           .... *
// *       </TaskSequence> *
// * *
// *    </KBFit> *
// * *
// * *
// *   (a) If <ProjectName> is missing, a default name will be created. *
// * *
// *   (b) If <Logfile> is missing, a default name for the log file is used. *
// * *
// *   (c) If <EchoXML> is missing, the input XML will not be written to the *
// *       log file. *
// * *
// *   (d) The tag <MCSamplingInfo> is mandatory.  It controls the default *
// *       resampling method:  jackknife or bootstrap.  This default method *
// *       is assumed for all reading and writing sampling results to and *
// *       from files.  Note that both jackknife and bootstrap resampling *
// *       can be done in any program execution, but only one can be used *
// *       for reading/writing to files.  This tag has the form below.  See *
// *       comments for the MCSamplingInfo and Bootstrapper classes for more *
// *       details about this tag. *
// * *
// *      <MCSamplingInfo> *
// *         <Jackknife/> *
// *      </MCSamplingInfo> *
// *                       OR *
// *      <MCSamplingInfo> *
// *         <Bootstrapper> *
// *            <NumberResamplings>2048</NumberResamplings> *
// *            <Seed>6754</Seed> *
// *            <BootSkip>127</BootSkip> *
// *         </Bootstrapper> *
// *      </MCSamplingInfo> *
// * *
// ******************************************************************************

int main(int argc, const char* argv[]) {
  // convert arguments to C++ strings
  vector<string> tokens(argc - 1);
  for (int k = 1; k < argc; ++k) {
    tokens[k - 1] = string(argv[k]);
  }

  if (tokens.size() != 1) {
    cout << "Error: requires a file name as the only argument" << endl;
    return 1;
  }

  try {
    XMLHandler xmltask;
    xmltask.set_exceptions_on();
    if (tokens.size() > 0) {
      string filename(tokens[0]);
      xmltask.set_from_file(filename);
    }

    // set up the task handler
    TaskHandler tasker(xmltask);

    // do the tasks in sequence
    tasker.do_batch_tasks(xmltask);
  } catch (const std::exception& msg) {
    cout << "Error: " << msg.what() << endl;
    return 1;
  }

  return 0;
}
