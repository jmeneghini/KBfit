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

// Function to display help information
void show_help() {
  cout << endl;
  cout << "=====================================================" << endl;
  cout << "                     KBfit                           " << endl;
  cout << "=====================================================" << endl;
  cout << endl;
  cout << "DESCRIPTION:" << endl;
  cout << "  KBfit is a lattice QCD analysis tool for studying two-hadron" << endl;
  cout << "  interactions in a finite box using the Luscher method." << endl;
  cout << "  It performs K-matrix fitting, along with helpers to assist in" << endl;
  cout << "  the analysis." << endl;
  cout << endl;
  cout << "USAGE:" << endl;
  cout << "  KBfit <input_file.xml>" << endl;
  cout << "  KBfit -h | --help" << endl;
  cout << endl;
  cout << "ARGUMENTS:" << endl;
  cout << "  input_file.xml    XML input file containing analysis configuration" << endl;
  cout << "  -h, --help        Show this help message and exit" << endl;
  cout << endl;
  cout << "AVAILABLE TASKS:" << endl;
  cout << "  DoPrint           Print quantization condition results and spectra" << endl;
  cout << "  DoFit             Perform K-matrix fitting to lattice data" << endl;
  cout << "  DoSingleChannel   Single-channel scattering analysis" << endl;
  cout << endl;
  cout << "XML INPUT STRUCTURE:" << endl;
  cout << "  The input XML file must have the following structure:" << endl;
  cout << endl;
  cout << "  <KBFit>" << endl;
  cout << "    <Initialize>" << endl;
  cout << "      <ProjectName>YourProjectName</ProjectName>" << endl;
  cout << "      <LogFile>output.log</LogFile>              [optional]" << endl;
  cout << "      <EchoXML/>                                 [optional]" << endl;
  cout << "      <MCSamplingInfo>                           [required]" << endl;
  cout << "        <Jackknife/> OR <Bootstrapper>...</Bootstrapper>" << endl;
  cout << "      </MCSamplingInfo>" << endl;
  cout << "    </Initialize>" << endl;
  cout << endl;
  cout << "    <TaskSequence>" << endl;
  cout << "      <Task>" << endl;
  cout << "        <Action>DoPrint|DoFit|DoSingleChannel</Action>" << endl;
  cout << "        <!-- Task-specific configuration -->" << endl;
  cout << "      </Task>" << endl;
  cout << "      <!-- Additional tasks as needed -->" << endl;
  cout << "    </TaskSequence>" << endl;
  cout << "  </KBFit>" << endl;
  cout << endl;
  cout << "NOTES:" << endl;
  cout << "  - If <ProjectName> is missing, a default name will be created" << endl;
  cout << "  - If <LogFile> is missing, a default log file name is used" << endl;
  cout << "  - <MCSamplingInfo> is mandatory and controls the resampling method" << endl;
  cout << "  - <EchoXML/> causes the input XML to be written to the log file" << endl;
  cout << "  - Example XML files are available in the examples/ directory" << endl;
  cout << endl;
  cout << "OUTPUT:" << endl;
  cout << "  - Log file with detailed analysis results" << endl;
  cout << "  - CSV files with numerical results (task-dependent)" << endl;
  cout << "  - HDF5 files for sampling results (when applicable)" << endl;
  cout << endl;
  cout << "For more detailed information, please refer to the" << endl;
  cout << "example XML files in the examples/ directory." << endl;
  cout << endl;
}

int main(int argc, const char* argv[]) {
  // convert arguments to C++ strings
  vector<string> tokens(argc - 1);
  for (int k = 1; k < argc; ++k) {
    tokens[k - 1] = string(argv[k]);
  }

  // Check for help arguments
  if (tokens.size() == 1 && (tokens[0] == "-h" || tokens[0] == "--help")) {
    show_help();
    return 0;
  }

  // Check for correct number of arguments
  if (tokens.size() != 1) {
    cout << "Error: KBfit requires exactly one argument (the XML input file)" << endl;
    cout << "Usage: KBfit <input_file.xml>" << endl;
    cout << "       KBfit -h | --help" << endl;
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
