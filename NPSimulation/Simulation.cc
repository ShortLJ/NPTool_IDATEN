//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// STL
#include<cstdlib>

#include "G4RunManager.hh"
#include "G4PhysListFactory.hh"
// UI
#include "G4UImanager.hh"
#include "QBBC.hh"

#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4VisManager.hh"

#include "G4VisExecutive.hh"

#include "G4UIExecutive.hh"

// G4 local source
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

// G4 General Source
#include "SteppingVerbose.hh"
#include "Randomize.hh"

// Root
#include "TRandom.h"

// NPS headers
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "NPSimulationVersion.hh"
//NPL headers
#include "NPOptionManager.h"
#include "RootOutput.h"

int main(int argc, char** argv){
    // Initialize NPOptionManager object
    NPOptionManager* OptionManager  = NPOptionManager::getInstance(argc, argv);
    OptionManager->SetIsSimulation();
    if(OptionManager->GetVerboseLevel() > 0){
        string line;
        line.resize(80,'*');
        cout << endl << line << endl;
        cout << "********************************  NPSimulation  ********************************"<< endl;
        cout << line << endl;
        cout << "NPSimulation version: npsimulation-"<< NPS::version_major <<"-" << NPS::version_minor << "-" << NPS::version_dets <<endl;
        cout << " Copyright: NPTool Collaboration "<<endl;
        cout << " Gitlab: https://gitlab.in2p3.fr/np/nptool "<<endl; ;
        cout << line << endl;
    }
    
    // Test if input files are found. If not, exit
    if (OptionManager->IsDefault("EventGenerator"))
        OptionManager->SendErrorAndExit("EventGenerator");
    if (OptionManager->IsDefault("DetectorConfiguration"))
        OptionManager->SendErrorAndExit("DetectorConfiguration");
    // case when input files are here
    G4String EventGeneratorFileName = OptionManager->GetReactionFile();
    G4String DetectorFileName       = OptionManager->GetDetectorFile();
    
    
    // initialize the state of the root and geant4 random generator
    if(OptionManager->GetRandomSeed()>0){
      std::cout << " Seeds for random generators set to: " << OptionManager->GetRandomSeed() << std::endl;
      gRandom->SetSeed(OptionManager->GetRandomSeed()); 
      CLHEP::HepRandom::setTheSeed(OptionManager->GetRandomSeed(),3);
    }
    
    
    
    // my Verbose output class
    G4VSteppingVerbose::SetInstance(new SteppingVerbose);
    
    ///////////////////////////////////////////////////////////////
    ///////////////// Initializing the Root Output ////////////////
    ///////////////////////////////////////////////////////////////
    RootOutput::getInstance(OptionManager->GetOutputFile());
    
    // Construct the default run manager
    G4RunManager* runManager = new G4RunManager;
    
    // set mandatory initialization classes
    DetectorConstruction* detector  = new DetectorConstruction();
    runManager->SetUserInitialization(detector);
    
    PhysicsList* physicsList   = new PhysicsList();
    physicsList->SetVerboseLevel(0);
    runManager->SetUserInitialization(physicsList);
    PrimaryGeneratorAction* primary = new PrimaryGeneratorAction(detector);

    // Initialize Geant4 kernel
    runManager->Initialize();
    
    ///////////////////////////////////////////////////////////////
    /////////// Define UI terminal for interactive mode ///////////
    ///////////////////////////////////////////////////////////////
    // interactive mode : define UI session
    // Get the pointer to the User Interface manager
    G4cout << "//////////// Starting UI ////////////"<< endl;
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    
    
    ///////////////////////////////////////////////////////////////
    ////////////////////// Reading Reaction ///////////////////////
    ///////////////////////////////////////////////////////////////
    primary->ReadEventGeneratorFile(EventGeneratorFileName);
    runManager->SetUserAction(primary);
  
    ///////////////////////////////////////////////////////////////
    ///////////////// Starting the Stepping Action ////////////////
    ///////////////////////////////////////////////////////////////
    SteppingAction* stepping_action = new SteppingAction() ;
    runManager->SetUserAction(stepping_action)       ;
    
    ///////////////////////////////////////////////////////////////
    ////////////////// Starting the Event Action //////////////////
    ///////////////////////////////////////////////////////////////
    EventAction* event_action = new EventAction() ;
    event_action->SetDetector(detector)           ;
    runManager->SetUserAction(event_action)       ;
    
    ///////////////////////////////////////////////////////////////
    /////////////////// Starting the Run Action ///////////////////
    ///////////////////////////////////////////////////////////////
    RunAction* run_action = new RunAction() ;
    runManager->SetUserAction(run_action);
    
    
    G4VisManager* visManager=NULL;

    if(!OptionManager->GetG4BatchMode()){
        string Path_Macro = getenv("NPTOOL");
        Path_Macro+="/NPSimulation/ressources/macro/";
        UImanager->ApplyCommand("/control/execute " +Path_Macro+"verbose.mac");

        UImanager->ApplyCommand("/control/execute " +Path_Macro+"aliases.mac");
        visManager = new G4VisExecutive("Quiet");
        visManager->Initialize();
        UImanager->ApplyCommand("/control/execute " +Path_Macro+"vis.mac");
        if (ui->IsGUI()){
            UImanager->ApplyCommand("/control/execute " +Path_Macro+"gui.mac");
        }

  #ifdef __APPLE__
        string command= "osascript ";
        command+= getenv("NPTOOL");
        command+="/NPSimulation/ressources/scripts/bringtofront.osa & ";
        int res =system(command.c_str());
        res =0;
  #endif
    }
    else{// if batch mode do not accumulate any track
        UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate 0");
      }
    // Execute user macro
    if(!OptionManager->IsDefault("G4MacroPath")){
        UImanager->ApplyCommand("/control/execute "+ OptionManager->GetG4MacroPath());
    }
    
    // Start the session
    if(!OptionManager->GetG4BatchMode())
        ui->SessionStart();
    
    
    delete ui;
    
    if(visManager)
        delete visManager;
    
    ///////////////////////////////////////////////////////////////
    ////////////////////// Job termination ////////////////////////
    ///////////////////////////////////////////////////////////////
    // Save the Geant4 random generator internal generator state in a TASCII 
    // file store with the root output
    std::ofstream file(".geant4_random_state");
    CLHEP::HepRandom::saveFullState(file);  
    file.close(); 
    TAsciiFile* aFile = new TAsciiFile();
    aFile->SetNameTitle("G4RandomFinalState",".geant4_random_state");
    aFile->Append(".geant4_random_state");
    RootOutput::getInstance()->GetFile()->cd();
    aFile->Write(0,TAsciiFile::kOverwrite);
    int dummy = system("rm .geant4_random_state");
    dummy*=2;
    // delete primary; delete detector;
    
    delete runManager;
    
    RootOutput::getInstance()->Destroy();
    return 0;
}
