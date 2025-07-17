#ifndef TASK_HANDLER_H
#define TASK_HANDLER_H

#include "kbobs_handler.h"
#include "xml_handler.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <mpi.h>

// ***************************************************************

//  Create a directory with the project name and quantization condition
//  if it does not already exist.

std::filesystem::path
createKBOutputDirectory(const std::string& base_directory = "",
                        const std::string& quantization_condition = "");

class TaskHandlerData; // base class for persistent data

// ******************************************************************************
// * *
// *   "TaskHandler" is the workhorse class of KBfit.  It manages the tasks *
// *   to do, and maintains pointers to objects that get the data and analyze *
// *   the data.  The constructor sets up the known tasks and creates the *
// *   data handling objects.  An object of class "KBObsHandler" is created *
// *   and is accessed through the pointer "m_obs"; this object is *
// *   responsible for reading data from files.
// * *
// *   The "TaskHandler" constructor requires XML of the form *
// * *
// *    <KBFit> *
// * *
// *       <Initialize> *
// *         <ProjectName>NameOfProject</ProjectName> *
// *         <OutputDirectory>output_directory</OutputDirectory> *
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
// *   (b) If <OutputDirectory> is missing, the current working directory *
// *       is used. *
// * *
// *   (c) If <LogFile> is missing, a default name for the log file is used. *
// *       The log file is put in the output directory if no path separators *
// *       are found in the log file name.                                   *
// * *
// *   (d) If <EchoXML> is missing, the input XML will not be written to the *
// *       log file. *
// * *
// *   (e) The tag <MCSamplingInfo> is mandatory.  It controls the default *
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
// * *
// ******************************************************************************

class TaskHandler {

  KBObsHandler* m_obs;
  std::ofstream clog;
  std::string m_project_name;
  std::string m_output_directory;
  int m_mpi_rank; // MPI rank for this process

  typedef void (TaskHandler::*task_ptr)(XMLHandler&, XMLHandler&, int);
  std::map<std::string, task_ptr> m_task_map;
  std::map<std::string, TaskHandlerData*> m_task_data_map;

  // Prevent copying ... handler might contain large
  // amounts of data

#ifndef NO_CXX11
  TaskHandler() = delete;
  TaskHandler(const TaskHandler&) = delete;
  TaskHandler& operator=(const TaskHandler&) = delete;
#else
  TaskHandler();
  TaskHandler(const TaskHandler&);
  TaskHandler& operator=(const TaskHandler&);
#endif

public:
  TaskHandler(XMLHandler& xmlin, int mpi_rank = 0);
  ~TaskHandler();

  void do_batch_tasks(XMLHandler& xmlin);

  void clear_task_data();
  void erase_task_data(const std::string& tdname);
  void insert_task_data(const std::string& tdname, TaskHandlerData* tdata);
  TaskHandlerData* get_task_data(const std::string& tdname);
  KBObsHandler* getKBObsHandler() { return m_obs; }

private:
  void do_task(XMLHandler& xml_in, XMLHandler& output, int taskcount);

  std::string get_date_time();

  void finish_log();

  // The important task subroutines

  void doPrint(XMLHandler& xml_in, XMLHandler& output, int taskcount);
  void doFit(XMLHandler& xml_in, XMLHandler& output, int taskcount);
  void doSingleChannel(XMLHandler& xml_in, XMLHandler& output, int taskcount);
};

// ***************************************************************

// base class for persistent data

class TaskHandlerData {
public:
  TaskHandlerData() {}
  virtual ~TaskHandlerData() {}
};

// ***************************************************************
#endif
