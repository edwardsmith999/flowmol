// ATC transfer headers
#include "TimeIntegrator.h"
#include "ATC_Coupling.h"
#include "ATC_Error.h"

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeIntegrator
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  TimeIntegrator::TimeIntegrator(ATC_Coupling * atc,
                                 TimeIntegrationType timeIntegrationType) :
    timeIntegrationMethod_(NULL),
    atc_(atc),
    timeFilter_(NULL),
    timeFilterManager_(atc_->time_filter_manager()),
    timeIntegrationType_(timeIntegrationType),
    needReset_(true)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  TimeIntegrator::~TimeIntegrator()
  {
    if (timeIntegrationMethod_)
      delete timeIntegrationMethod_;
  }

  //--------------------------------------------------------
  //  construct_transfers
  //    sets/constructs all required dependency managed data
  //--------------------------------------------------------
  void TimeIntegrator::construct_transfers()
  {
    timeIntegrationMethod_->construct_transfers();
  }

  //--------------------------------------------------------
  //  initialize
  //    initialize all data and variables before a run
  //--------------------------------------------------------
  void TimeIntegrator::initialize()
  {
    timeIntegrationMethod_->initialize();
    needReset_ = false;
  }

  //--------------------------------------------------------
  //  pre_initial_integrate1
  //    first time integration computations
  //    before Verlet step 1
  //--------------------------------------------------------
  void TimeIntegrator::pre_initial_integrate1(double dt)
  {
    timeIntegrationMethod_->pre_initial_integrate1(dt);
  }

  //--------------------------------------------------------
  //  pre_initial_integrate2
  //    second time integration computations
  //    before Verlet step 1
  //--------------------------------------------------------
  void TimeIntegrator::pre_initial_integrate2(double dt)
  {
    timeIntegrationMethod_->pre_initial_integrate2(dt);
  }

  //--------------------------------------------------------
  //  mid_initial_integrate1
  //    first time integration computations
  //    at the mid-point of Verlet step 1
  //--------------------------------------------------------
  void TimeIntegrator::mid_initial_integrate1(double dt)
  {
    timeIntegrationMethod_->mid_initial_integrate1(dt);
  }

  //--------------------------------------------------------
  //  mid_initial_integrate2
  //    second time integration computations
  //    at the mid-point of Verlet step 1
  //--------------------------------------------------------
  void TimeIntegrator::mid_initial_integrate2(double dt)
  {
    timeIntegrationMethod_->mid_initial_integrate2(dt);
  }

  //--------------------------------------------------------
  //  post_initial_integrate1
  //    first time integration computations
  //    after Verlet step 1
  //--------------------------------------------------------
  void TimeIntegrator::post_initial_integrate1(double dt)
  {
    timeIntegrationMethod_->post_initial_integrate1(dt);
  }

  //--------------------------------------------------------
  //  post_initial_integrate2
  //    second time integration computations
  //    after Verlet step 1
  //--------------------------------------------------------
  void TimeIntegrator::post_initial_integrate2(double dt)
  {
    timeIntegrationMethod_->post_initial_integrate2(dt);
  }

  //--------------------------------------------------------
  //  pre_final_integrate1
  //    first time integration computations 
  //    before Verlet step 2
  //--------------------------------------------------------
  void TimeIntegrator::pre_final_integrate1(double dt)
  {
    timeIntegrationMethod_->pre_final_integrate1(dt);
  }

  //--------------------------------------------------------
  //  pre_final_integrate2
  //    second time integration computations
  //    before Verlet step 2
  //--------------------------------------------------------
  void TimeIntegrator::pre_final_integrate2(double dt)
  {
    timeIntegrationMethod_->pre_final_integrate2(dt);
  }

  //--------------------------------------------------------
  //  post_final_integrate1
  //    first time integration computations 
  //    after Verlet step 2
  //--------------------------------------------------------
  void TimeIntegrator::post_final_integrate1(double dt)
  {
    timeIntegrationMethod_->post_final_integrate1(dt);
  }

  //--------------------------------------------------------
  //  post_final_integrate2
  //    second time integration computations
  //    after Verlet step 2
  //--------------------------------------------------------
  void TimeIntegrator::post_final_integrate2(double dt)
  {
    timeIntegrationMethod_->post_final_integrate2(dt);
  }

  //--------------------------------------------------------
  //  post_final_integrate3
  //    third time integration computations
  //    after Verlet step 2
  //--------------------------------------------------------
  void TimeIntegrator::post_final_integrate3(double dt)
  {
    timeIntegrationMethod_->post_final_integrate3(dt);
  }

  //--------------------------------------------------------
  //  has_final_predictor
  //    checks to see if first RHS computation is needed
  //--------------------------------------------------------
  bool TimeIntegrator::has_final_predictor()
  {
    return timeIntegrationMethod_->has_final_predictor();
  }

  //--------------------------------------------------------
  //  has_final_corrector
  //    checks to see if second RHS computation is needed
  //--------------------------------------------------------
  bool TimeIntegrator::has_final_corrector()
  {
    return timeIntegrationMethod_->has_final_corrector();
  }

  //--------------------------------------------------------
  //  add_to_rhs
  //    add any needed contributions to RHS
  //--------------------------------------------------------
  void TimeIntegrator::add_to_rhs()
  {
    timeIntegrationMethod_->add_to_rhs();
  }

  //--------------------------------------------------------
  //  post_process
  //    perform any post processing calculations
  //--------------------------------------------------------
  void TimeIntegrator::post_process()
  {
    timeIntegrationMethod_->post_process();
  }

  //--------------------------------------------------------
  //  output
  //    add variables to output list
  //--------------------------------------------------------
  void TimeIntegrator::output(OUTPUT_LIST & outputData)
  {
    timeIntegrationMethod_->output(outputData);
  }

  //--------------------------------------------------------
  //  pack_fields
  //    add persistent variables to data list
  //--------------------------------------------------------
  void TimeIntegrator::pack_fields(RESTART_LIST & data)
  {
    timeIntegrationMethod_->pack_fields(data);
    
    //timeFilter_->pack_fields(data);
  }

  //--------------------------------------------------------
  //  finish
  //    perform any final state setting
  //--------------------------------------------------------
  void TimeIntegrator::finish()
  {
    timeIntegrationMethod_->finish();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeIntegrationMethod
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //        Grab data from ATC
  //-------------------------------------------------------- 
  TimeIntegrationMethod::TimeIntegrationMethod(TimeIntegrator * timeIntegrator) :
    timeIntegrator_(timeIntegrator),
    atc_(timeIntegrator_->atc())
  {
    // do nothing
  }

};
