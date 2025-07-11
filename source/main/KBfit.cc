#include "task_handler.h"
#include <cstdio>
#include <ctime>
#include <iostream>
#include <map>
#include <vector>
#include <mpi.h>

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
  // Check for MPI worker flag before MPI initialization
  bool is_mpi_worker = false;
  for (int i = 1; i < argc; ++i) {
    if (string(argv[i]) == "--mpi-worker") {
      is_mpi_worker = true;
      break;
    }
  }
  
  // Initialize MPI
  MPI_Init(&argc, const_cast<char***>(&argv));
  
  // Handle MPI worker processes
  if (is_mpi_worker) {
    // This is a spawned MPI worker process
    std::cerr << "MPI worker process started" << std::endl;
    
    // The worker process will participate in the merged communicator
    // and then exit. The actual work is handled in doChiSquareFittingMPI
    // which is called through the MPI communication protocol.
    
    // For now, just participate in the communicator and exit
    MPI_Comm parent_comm;
    MPI_Comm_get_parent(&parent_comm);
    
    if (parent_comm != MPI_COMM_NULL) {
      // Merge with parent communicator
      MPI_Comm merged_comm;
      MPI_Intercomm_merge(parent_comm, 1, &merged_comm);  // 1 = high group (child)
      
      int rank, size;
      MPI_Comm_rank(merged_comm, &rank);
      MPI_Comm_size(merged_comm, &size);
      
      std::cerr << "MPI worker process " << rank << " of " << size << " ready" << std::endl;
      
      // Receive basic parameters to participate in broadcasts
      uint params_data[4];
      MPI_Bcast(params_data, 4, MPI_UNSIGNED, 0, merged_comm);
      
      // Receive other broadcast data
      double chisq_dof_recv;
      MPI_Bcast(&chisq_dof_recv, 1, MPI_DOUBLE, 0, merged_comm);
      
      uint nparams = params_data[0];
      vector<double> params_fullsample(nparams);
      MPI_Bcast(params_fullsample.data(), nparams, MPI_DOUBLE, 0, merged_comm);
      
      // Receive minimizer info
      int csm_str_len;
      MPI_Bcast(&csm_str_len, 1, MPI_INT, 0, merged_comm);
      string csm_str(csm_str_len, '\0');
      MPI_Bcast(const_cast<char*>(csm_str.c_str()), csm_str_len, MPI_CHAR, 0, merged_comm);
      
      std::cerr << "MPI worker " << rank << " received parameters (serialization not implemented)" << std::endl;
      
      // Clean up and exit
      MPI_Comm_free(&merged_comm);
      std::cerr << "MPI worker " << rank << " finished" << std::endl;
    } else {
      std::cerr << "Error: MPI worker process has no parent communicator" << std::endl;
    }
    
    MPI_Finalize();
    return 0;
  }
  
  // Get MPI info for logging (should be size 1 for main process)
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  // Main program should run as single process
  if (size > 1) {
    if (rank == 0) {
      cout << "Warning: KBfit detected " << size << " MPI processes at startup." << endl;
      cout << "Only rank 0 will run the main program. MPI processes will be spawned" << endl;
      cout << "dynamically for chi-square fitting when needed." << endl;
    } else {
      // Non-zero ranks should exit immediately
      MPI_Finalize();
      return 0;
    }
  }
  
  // Only rank 0 continues with the main program
  rank = 0;  // Force rank 0 for the main program logic
  
  // Only rank 0 handles command line processing and output
  int return_code = 0;
  
  cout << "KBfit starting with dynamic MPI process spawning for parallel fitting" << endl;
  
  // convert arguments to C++ strings
  vector<string> tokens(argc - 1);
  for (int k = 1; k < argc; ++k) {
    tokens[k - 1] = string(argv[k]);
  }

  // Check for help arguments
  if (tokens.size() == 1 && (tokens[0] == "-h" || tokens[0] == "--help")) {
    show_help();
    MPI_Finalize();
    return 1;
  }
  // Check for correct number of arguments
  else if (tokens.size() != 1) {
    cout << "Error: KBfit requires exactly one argument (the XML input file)" << endl;
    cout << "Usage: KBfit <input_file.xml>" << endl;
    cout << "       KBfit -h | --help" << endl;
    MPI_Finalize();
    return 1;
  }

  // Run the main program (single process)
  try {
    XMLHandler xmltask;
    xmltask.set_exceptions_on();
    
    // Read the input file
    string filename(tokens[0]);
    xmltask.set_from_file(filename);
    
    // Set up the task handler (only rank 0)
    TaskHandler tasker(xmltask, 0);

    // Do the tasks in sequence
    tasker.do_batch_tasks(xmltask);
    return_code = 0;
  } catch (const std::exception& msg) {
    cout << "Error: " << msg.what() << endl;
    return_code = 1;
  }

  // Finalize MPI
  MPI_Finalize();
  
  return return_code;
}
