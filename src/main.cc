/*  Copyright 2011-2013 Alexis Herault, Giuseppe Bilotta, Robert A. Dalrymple, Eugenio Rustico, Ciro Del Negro

    Istituto Nazionale di Geofisica e Vulcanologia
        Sezione di Catania, Catania, Italy

    Università di Catania, Catania, Italy

    Johns Hopkins University, Baltimore, MD

    This file is part of GPUSPH.

    GPUSPH is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GPUSPH is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GPUSPH.  If not, see <http://www.gnu.org/licenses/>.
*/


// endl, cerr, etc.
#include <iostream>
// signal, sigaction, etc.
#include <signal.h>

#include "GPUSPH.h"
#include "Options.h"
#include "GlobalData.h"
#include "NetworkManager.h"

// Include only the problem selected at compile time (PROBLEM, QUOTED_PROBLEM)
#include "problem_select.opt"

/* Include all other opt file for show_version */
#include "gpusph_version.opt"
#include "fastmath_select.opt"
#include "dbg_select.opt"
#include "compute_select.opt"

using std::cout;
using std::endl;
using std::invalid_argument;
using std::cerr;

void show_version()
{
	static const char dbg_or_rel[] =
#if defined(_DEBUG_)
		"Debug";
#else
		"Release";
#endif

	printf("GPUSPH version %s\n", GPUSPH_VERSION);
	printf("%s version %s fastmath for compute capability %u.%u\n",
		dbg_or_rel,
		FASTMATH ? "with" : "without",
		COMPUTE/10, COMPUTE%10);
	printf("Compiled for problem \"%s\"\n", QUOTED_PROBLEM);
}


// TODO: cleanup, no exit
void print_usage() {
	show_version();
	cout << "Syntax: " << endl;
	cout << "\tGPUSPH [--device n[,n...]] [--dem dem_file] [--deltap VAL] [--tend VAL]\n";
	cout << "\t       [--resume fname] [--checkpoint-every VAL] [--checkpoints VAL]\n";
	cout << "\t       [--dir directory] [--nosave] [--striping] [--gpudirect [--asyncmpi]]\n";
	cout << "\t       [--num-hosts VAL [--byslot-scheduling]]\n";
	cout << "\tGPUSPH --help\n\n";
	cout << " --resume : resume from the given file (HotStart file saved by HotWriter)\n";
	cout << " --checkpoint-every : HotStart checkpoints will be created every VAL seconds\n";
	cout << "                      of simulated time (float VAL, 0 disables)\n";
	cout << " --checkpoints : number of HotStart checkpoints to keep (integer VAL)\n";
	cout << " --device n[,n...] : Use device number n; runs multi-gpu if multiple n are given\n";
	cout << " --dem : Use given DEM (if problem supports it)\n";
	cout << " --deltap : Use given deltap (VAL is cast to float)\n";
	cout << " --tend : Break at given time (VAL is cast to float)\n";
	cout << " --maxiter : Break after this many iterations (integer VAL)\n";
	cout << " --dir : Use given directory for dumps instead of date-based one\n";
	cout << " --nosave : Disable all file dumps but the last\n";
	cout << " --gpudirect: Enable GPUDirect for RDMA (requires a CUDA-aware MPI library)\n";
	cout << " --striping : Enable computation/transfer overlap  in multi-GPU (usually convenient for 3+ devices)\n";
	cout << " --asyncmpi : Enable asynchronous network transfers (requires GPUDirect and 1 process per device)\n";
	cout << " --num-hosts : Uses multiple processes per node by specifying the number of nodes (VAL is cast to uint)\n";
	cout << " --byslot-scheduling : MPI scheduler is filling hosts first, as opposite to round robin scheduling\n";
	//cout << " --nobalance : Disable dynamic load balancing\n";
	//cout << " --lb-threshold : Set custom LB activation threshold (VAL is cast to float)\n";
	cout << " --help: Show this help and exit\n";
}

// if some option needs to be passed to GlobalData, remember to set it in GPUSPH::initialize()
int parse_options(int argc, char **argv, GlobalData *gdata)
{
	const char *arg(NULL);

	if (!gdata) return -1;
	Options* _clOptions = gdata->clOptions;

	// skip arg 0 (program name)
	argv++; argc--;

	while (argc > 0) {
		arg = *argv;
		argv++;
		argc--;
		if (!strcmp(arg, "--resume")) {
			_clOptions->resume_fname = std::string(*argv);
			argv++;
			argc--;
		} else if (!strcmp(arg, "--checkpoint-every")) {
			sscanf(*argv, "%f", &(_clOptions->checkpoint_freq));
			argv++;
			argc--;
		} else if (!strcmp(arg, "--checkpoints")) {
			sscanf(*argv, "%d", &(_clOptions->checkpoints));
			argv++;
			argc--;
		} else if (!strcmp(arg, "--device")) {
			/* read the next arg as a list of integers */
			char * pch;
			pch = strtok (*argv,",");
			while (pch != NULL) {
				//printf ("%s\n",pch);
				if (gdata->devices==MAX_DEVICES_PER_NODE) {
					printf("WARNING: devices exceeding number %u will be ignored\n",
						gdata->device[MAX_DEVICES_PER_NODE-1]);
					break;
				} else {
					// inc _clOptions->devices only if scanf was successful
					if (sscanf(pch, "%u", &(gdata->device[gdata->devices]))>0) {
						gdata->devices++;
						gdata->totDevices++;
					} else {
						printf("WARNING: token %s is not a number - ignored\n", pch);
						//break;
					}
				}
				pch = strtok (NULL, " ,.-");
			}
			if (gdata->devices<1) {
				fprintf(stderr, "ERROR: --device option given, but no device specified\n");
				return -1;
			}
			argv++;
			argc--;
		} else if (!strcmp(arg, "--deltap")) {
			/* read the next arg as a double */
			sscanf(*argv, "%lf", &(_clOptions->deltap));
			argv++;
			argc--;
		} else if (!strcmp(arg, "--tend")) {
			/* read the next arg as a float */
			sscanf(*argv, "%f", &(_clOptions->tend));
			argv++;
			argc--;
		} else if (!strcmp(arg, "--maxiter")) {
			/* read the next arg as a int */
			sscanf(*argv, "%d", &(_clOptions->maxiter));
			argv++;
			argc--;
		} else if (!strcmp(arg, "--dem")) {
			_clOptions->dem = std::string(*argv);
			argv++;
			argc--;
		} else if (!strcmp(arg, "--dir")) {
			_clOptions->dir = std::string(*argv);
			argv++;
			argc--;
		} else if (!strcmp(arg, "--nosave")) {
			_clOptions->nosave = true;
		} else if (!strcmp(arg, "--gpudirect")) {
			_clOptions->gpudirect = true;
		} else if (!strcmp(arg, "--striping")) {
			_clOptions->striping = true;
		} else if (!strcmp(arg, "--asyncmpi")) {
			_clOptions->asyncNetworkTransfers = true;
		} else if (!strcmp(arg, "--num-hosts") || !strcmp(arg, "--num_hosts")) {
			/* read the next arg as a uint */
			sscanf(*argv, "%u", &(_clOptions->num_hosts));
			argv++;
			argc--;
		} else if (!strcmp(arg, "--version")) {
			show_version();
			return 0;
		} else if (!strcmp(arg, "--byslot-scheduling") || !strcmp(arg, "--byslot_scheduling")) {
			_clOptions->byslot_scheduling = true;
#if 0 // options will be enabled later
		} else if (!strcmp(arg, "--nobalance")) {
			_clOptions->nobalance = true;
		} else if (!strcmp(arg, "--lb-threshold")) {
			// read the next arg as a float
			sscanf(*argv, "%f", &(_clOptions->custom_lb_threshold));
			argv++;
			argc--;
#endif
		} else if (!strcmp(arg, "--help")) {
			print_usage();
			return 0;
		} else if (!strncmp(arg, "--", 2)) {
			// TODO bool options would need to be treated specially,
			// currently they require a following 1 or 0
			_clOptions->set(arg+2, *argv);
			argv++;
			argc--;
		} else {
			cerr << "Fatal: Unknown option: " << arg << endl;
			return -1;

			// Left for future dynamic loading:
			/*if (_clOptions->problem.empty()) {
				_clOptions->problem = std::string(arg);
			} else {
				cout << "Problem " << arg << " selected after problem " << _clOptions->problem << endl;
			}*/
		}
	}

	if (gdata->devices==0) {
		printf(" * No devices specified, falling back to default (dev 0)...\n");
		// default: use first device. May use cutGetMaxGflopsDeviceId() instead.
		gdata->device[gdata->devices++] = 0;
	}

	// only for single-gpu
	_clOptions->device = gdata->device[0];

	_clOptions->problem = std::string( QUOTED_PROBLEM );

	// Left for future dynamic loading:
	/*if (_clOptions->problem.empty()) {
		problem_list();
		exit(0);
	}*/

	return 1;
}

bool check_short_length() {
	return (sizeof(uint) == 2*sizeof(short));
}

// SIGINT handler: issues a quit_request
void sigint_handler(int signum) {
	// if issued the second time, brutally terminate everything
	if (gdata_static_pointer->quit_request) {
		uint reachedt = gdata_static_pointer->threadSynchronizer->queryReachedThreads();
		uint maxt = gdata_static_pointer->threadSynchronizer->getNumThreads();
		if (reachedt > 0 && reachedt < maxt && !gdata_static_pointer->threadSynchronizer->didForceUnlockOccurr()) {
			printf("Second quit request - threads waiting: %u/%u. Forcing unlock...\n", reachedt, maxt);
			gdata_static_pointer->threadSynchronizer->forceUnlock();
		} else {
			printf("Unable to force unlock. Issuing exit(1)\n");
			exit(1);
		}
	} else {
		printf("Kindly asking to quit - please wait for the current iteration to complete\n");
		gdata_static_pointer->quit_request = true;
	}
}

void sigusr1_handler(int signum) {
	gdata_static_pointer->save_request = true;
}

int main(int argc, char** argv) {
	if (!check_short_length()) {
		printf("Fatal: this architecture does not have uint = 2 short\n");
		exit(1);
	}

	GlobalData gdata;
	gdata_static_pointer = &gdata;

	// Command line options
	gdata.clOptions = new Options();

	// catch SIGINT and SIGUSR1
	struct sigaction int_action, usr1_action;

	memset(&int_action, 0, sizeof(struct sigaction));
	int_action.sa_handler = sigint_handler;
	sigaction(SIGINT, &int_action, NULL);

	memset(&usr1_action, 0, sizeof(struct sigaction));
	usr1_action.sa_handler = sigusr1_handler;
	sigaction(SIGUSR1, &usr1_action, NULL);

	// parse command-line options
	int ret = parse_options(argc, argv, &gdata);
	if (ret <= 0)
		exit(ret);

	show_version();

	// TODO: check options, i.e. consistency

	// NOTE: Although GPUSPH has been designed to be run with one multi-threaded process per node, it is important not to create
	// any file or lock singleton resources before initializing the network, as the process might be forked
	gdata.networkManager = new NetworkManager();
	gdata.networkManager->initNetwork();
	gdata.networkManager->printInfo();

	int nm_worldsize = gdata.networkManager->getWorldSize();
	if (nm_worldsize > MAX_NODES_PER_CLUSTER) {
		cerr << "Too many nodes in cluster: " << nm_worldsize << " > " << MAX_NODES_PER_CLUSTER << endl;
		exit(1);
	}

	gdata.mpi_nodes = devcount_t(nm_worldsize);
	gdata.mpi_rank = gdata.networkManager->getProcessRank();

	// We "shift" the cuda device indices by devIndexOffset. It is useful in case of multiple processes per node. Will write external docs about the formula
	uint devIndexOffset = 0;
	if (gdata.clOptions->num_hosts > 0) {
		if (gdata.clOptions->byslot_scheduling)
			// non round-robin scheduling: fill first node, then start assigning to the second
			devIndexOffset = (gdata.mpi_rank % ( gdata.mpi_nodes / gdata.clOptions->num_hosts ) ) * gdata.devices;
		else
			// round-robin scheduling: distribute to non-empty node only if others have at least n-1 processes already
			devIndexOffset = (gdata.mpi_rank / gdata.clOptions->num_hosts) * gdata.devices;

		for (uint d=0; d < gdata.devices; d++)
				gdata.device[d] += devIndexOffset;
	} else
		if (gdata.clOptions->byslot_scheduling)
			printf("WARNING: --byslot-scheduling was enabled, but number of hosts is zero!\n");

	gdata.totDevices = gdata.mpi_nodes * gdata.devices;
	printf(" tot devs = %u (%u * %u)\n",gdata.totDevices, gdata.mpi_nodes, gdata.devices );
	if (gdata.clOptions->num_hosts > 0)
		printf(" num-hosts was specified: %u; shifting device numbers with offset %u\n", gdata.clOptions->num_hosts, devIndexOffset);

	if (gdata.clOptions->asyncNetworkTransfers) {

		if (!gdata.clOptions->gpudirect) {
			// since H2D and D2H transfers have to wait for network transfers
			fprintf(stderr, "FATAL: asynchronous network transfers require --gpudirect\n");
			gdata.networkManager->finalizeNetwork();
			return 1;
		}

		if (gdata.devices > 1) {
			// since we were too lazy to implement a more complex mechanism
			fprintf(stderr, "FATAL: asynchronous network transfers only supported with 1 process per device\n");
			gdata.networkManager->finalizeNetwork();
			return 1;
		}

	}

	// the Problem could (should?) be initialized inside GPUSPH::initialize()
	gdata.problem = new PROBLEM(&gdata);
	if (gdata.problem->m_simframework)
		gdata.simframework = gdata.problem->m_simframework;
	else
		throw invalid_argument("no simulation framework defined in the problem!");
	gdata.allocPolicy = gdata.simframework->getAllocPolicy();


	// get - and actually instantiate - the existing instance of GPUSPH
	GPUSPH *Simulator = GPUSPH::getInstance();

	// initialize CUDA, start workers, allocate CPU and GPU buffers
	bool initialized  = Simulator->initialize(&gdata);

	if (initialized) {
		printf("GPUSPH: initialized\n");

		// run the simulation until a quit request is triggered or an exception is thrown (TODO)
		Simulator->runSimulation();

		// finalize everything
		Simulator->finalize();
	} else
		printf("GPUSPH: problem during initialization, aborting...\n");

	// same consideration as above
	delete gdata.problem;

	// finalize MPI
	gdata.networkManager->finalizeNetwork();

	delete gdata.networkManager;

	return 0;
}

