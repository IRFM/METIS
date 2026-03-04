List of IMAS data produced by METIS
===================================


Introduction
------------

This document list all the possible IDS fields METIS can fill.
Some data are produced only when:
	-two point model is activated (egde data).
	-the model for tunsten is ativated.
	-the divertor contains tin.
	-plasma contains D & T.
	-plasma contains D & He3.
	-plasma contains H & B11.

The IDS equilibrium provides two 2D grids: one in cartesian coordinates (R,Z) ala EFIT and one inverse grid in (rho,theta).
IDSs named *_1 are alternative IDS for edge data stored in occurrence 1.
All other IDSs are stored in occurrence 0.

IDS pulse_schedule
-------------------

* pulse_schedule/time: Generic time
* pulse_schedule/ids_properties/comment: Any comment describing the content of this IDS
* pulse_schedule/ids_properties/homogeneous_time: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* pulse_schedule/ids_properties/source: Source of the data (any comment describing the origin of the data : code, path to diagnostic signals, processing method, ...). Superseeded by the new provenance structure.
* pulse_schedule/ids_properties/provider: Name of the person in charge of producing this data
* pulse_schedule/ids_properties/creation_date: Date at which this data has been produced
* pulse_schedule/density_control/mode: Control mode (operation mode and/or settings used by the controller)
* pulse_schedule/flux_control/mode: Control mode (operation mode and/or settings used by the controller)
* pulse_schedule/code/name: Name of software generating IDS
* pulse_schedule/code/commit: Unique commit reference of software
* pulse_schedule/code/version: Unique version (tag) of software
* pulse_schedule/code/repository: URL of software repository
* pulse_schedule/code/parameters: List of the code specific parameters in XML format
* pulse_schedule/code/output_flag: Output flag : 0 means the run is successful, other values mean some difficulty has been encountered, the exact meaning is then code specific. Negative values mean the result shall not be used.
* pulse_schedule/ic/power/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/ec/power_launched/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/lh/power/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/nbi/power/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/density_control/n_e_line/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/density_control/n_e_line/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/density_control/n_e_line/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/density_control/zeff/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/density_control/zeff/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/density_control/zeff/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/density_control/n_t_over_n_d/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/density_control/n_t_over_n_d/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/density_control/n_t_over_n_d/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/density_control/n_h_over_n_d/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/density_control/n_h_over_n_d/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/density_control/n_h_over_n_d/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/flux_control/i_plasma/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/flux_control/i_plasma/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/flux_control/i_plasma/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/flux_control/loop_voltage/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/flux_control/loop_voltage/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/flux_control/loop_voltage/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/flux_control/li_3/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/flux_control/li_3/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/flux_control/li_3/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/flux_control/beta_normal/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/flux_control/beta_normal/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/flux_control/beta_normal/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/position_control/minor_radius/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/minor_radius/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/minor_radius/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/position_control/elongation/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/elongation/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/elongation/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/position_control/elongation_upper/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/elongation_upper/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/elongation_upper/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/position_control/elongation_lower/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/elongation_lower/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/elongation_lower/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/position_control/triangularity/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/triangularity/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/triangularity/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/position_control/triangularity_upper/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/triangularity_upper/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/triangularity_upper/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/position_control/triangularity_lower/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/triangularity_lower/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/triangularity_lower/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/tf/b_field_tor_vacuum_r/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/tf/b_field_tor_vacuum_r/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/tf/b_field_tor_vacuum_r/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/ids_properties/provenance/node{:}/path: Path of the node within the IDS, following the syntax given in the link below. If empty, means the provenance information applies to the whole IDS.
* pulse_schedule/ic/antenna{:}/power/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/ic/antenna{:}/power/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/ic/antenna{:}/power/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/ic/antenna{:}/frequency/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/ic/antenna{:}/frequency/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/ic/antenna{:}/frequency/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/lh/antenna{:}/power/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/lh/antenna{:}/power/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/lh/antenna{:}/power/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/lh/antenna{:}/n_parallel/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/lh/antenna{:}/n_parallel/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/lh/antenna{:}/n_parallel/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/lh/antenna{:}/frequency/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/lh/antenna{:}/frequency/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/lh/antenna{:}/frequency/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/nbi/unit{:}/species{:}/label: String identifying neutral (e.g. H, D, T, He, C, ...)
* pulse_schedule/nbi/unit{:}/species{:}/name: String identifying neutral (e.g. H, D, T, He, C, ...)
* pulse_schedule/nbi/unit{:}/species{:}/fraction: Relative fraction of this species (in molecules) in the gas mixture
* pulse_schedule/nbi/unit{:}/power/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/nbi/unit{:}/power/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/nbi/unit{:}/power/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/nbi/unit{:}/energy/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/nbi/unit{:}/energy/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/nbi/unit{:}/energy/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/density_control/valve{:}/flow_rate/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/density_control/valve{:}/flow_rate/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/density_control/valve{:}/flow_rate/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/density_control/valve{:}/species{:}/label: String identifying neutral (e.g. H, D, T, He, C, ...)
* pulse_schedule/density_control/valve{:}/species{:}/name: String identifying neutral (e.g. H, D, T, He, C, ...)
* pulse_schedule/density_control/valve{:}/species{:}/fraction: Relative fraction of this species (in molecules) in the gas mixture
* pulse_schedule/position_control/magnetic_axis/r/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/magnetic_axis/r/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/magnetic_axis/r/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/position_control/magnetic_axis/z/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/magnetic_axis/z/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/magnetic_axis/z/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/position_control/geometric_axis/r/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/geometric_axis/r/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/geometric_axis/r/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/position_control/geometric_axis/z/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/geometric_axis/z/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/geometric_axis/z/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/position_control/x_point{:}/r/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/x_point{:}/r/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/x_point{:}/r/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/position_control/x_point{:}/z/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/x_point{:}/z/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/x_point{:}/z/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/position_control/boundary_outline{:}/r/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/boundary_outline{:}/r/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/boundary_outline{:}/r/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/position_control/boundary_outline{:}/z/reference_name: Reference name (e.g. in the native pulse schedule system of the device)
* pulse_schedule/position_control/boundary_outline{:}/z/reference: Reference waveform. Caution : error bars of the reference/data node are not used in the usual sense, instead they are used to describe the control envelope, with a meaning depending on the chosen envelope_type option.
* pulse_schedule/position_control/boundary_outline{:}/z/reference_type: Reference type:  0:relative (don't use for the moment, to be defined later when segments are introduced in the IDS structure); 1: absolute: the reference time trace is provided in the reference/data node
* pulse_schedule/ids_properties/provenance/node{:}/reference{:}/name: Reference name
* pulse_schedule/ids_properties/provenance/node{:}/reference{:}/timestamp: Date and time (UTC) at which the reference was created, expressed in a human readable form (ISO 8601) : the format of the string shall be : YYYY-MM-DDTHH:MM:SSZ. Example : 2020-07-24T14:19:00Z
* pulse_schedule/nbi/unit{:}/species{:}/element{:}/a: Mass of atom
* pulse_schedule/nbi/unit{:}/species{:}/element{:}/z_n: Nuclear charge
* pulse_schedule/density_control/valve{:}/species{:}/element{:}/a: Mass of atom
* pulse_schedule/density_control/valve{:}/species{:}/element{:}/z_n: Nuclear charge


IDS dataset_description
------------------------

* dataset_description/imas_version: Version of the IMAS infrastructure used to produce this data entry. Refers to the global IMAS repository which links to versions of every infrastructure tools
* dataset_description/dd_version: Version of the physics data dictionary of this dataset
* dataset_description/time: Generic time
* dataset_description/ids_properties/comment: Any comment describing the content of this IDS
* dataset_description/ids_properties/homogeneous_time: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* dataset_description/ids_properties/source: Source of the data (any comment describing the origin of the data : code, path to diagnostic signals, processing method, ...). Superseeded by the new provenance structure.
* dataset_description/ids_properties/provider: Name of the person in charge of producing this data
* dataset_description/ids_properties/creation_date: Date at which this data has been produced
* dataset_description/data_entry/user: Username
* dataset_description/data_entry/machine: Name of the experimental device to which this data is related
* dataset_description/data_entry/pulse_type: Type of the data entry, e.g. "pulse", "simulation", ...
* dataset_description/data_entry/pulse: Pulse number
* dataset_description/data_entry/run: Run number
* dataset_description/simulation/comment_before: Comment made when launching a simulation
* dataset_description/simulation/comment_after: Comment made at the end of a simulation
* dataset_description/simulation/time_begin: Start time
* dataset_description/simulation/time_step: Time interval between main steps, e.g. storage step (if relevant and constant)
* dataset_description/simulation/time_end: Stop time
* dataset_description/simulation/time_restart: Time of the last restart done during the simulation
* dataset_description/simulation/time_current: Current time of the simulation
* dataset_description/simulation/workflow: Description of the workflow which has been used to produce this data entry (e.g. copy of the Kepler MOML if using Kepler)
* dataset_description/ids_properties/provenance/node{:}/path: Path of the node within the IDS, following the syntax given in the link below. If empty, means the provenance information applies to the whole IDS.
* dataset_description/ids_properties/provenance/node{:}/reference{:}/name: Reference name
* dataset_description/ids_properties/provenance/node{:}/reference{:}/timestamp: Date and time (UTC) at which the reference was created, expressed in a human readable form (ISO 8601) : the format of the string shall be : YYYY-MM-DDTHH:MM:SSZ. Example : 2020-07-24T14:19:00Z


IDS dataset_fair
-----------------

* dataset_fair/identifier: Persistent identifier allowing to cite this data in a public and persistent way, should be provided as HTTP URIs
* dataset_fair/valid: Date range during which the data is or was valid. Expressed as YYYY-MM-DD/YYYY-MM-DD, where the former (resp. latter) date is the data at which the data started (resp. ceased) to be valid. If the data is still valid, the slash should still be present, i.e. indicate the validity start date with YYYY-MM-DD/. If the data ceased being valid but there is no information on the validity start date, indicate /YYYY-MM-DD.
* dataset_fair/rights_holder: The organisation owning or managing rights over this data
* dataset_fair/time: Generic time
* dataset_fair/ids_properties/comment: Any comment describing the content of this IDS
* dataset_fair/ids_properties/homogeneous_time: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* dataset_fair/ids_properties/source: Source of the data (any comment describing the origin of the data : code, path to diagnostic signals, processing method, ...). Superseeded by the new provenance structure.
* dataset_fair/ids_properties/provider: Name of the person in charge of producing this data
* dataset_fair/ids_properties/creation_date: Date at which this data has been produced
* dataset_fair/ids_properties/version_put/data_dictionary: Version of Data Dictionary used to PUT this IDS
* dataset_fair/ids_properties/version_put/access_layer: Version of Access Layer used to PUT this IDS
* dataset_fair/ids_properties/version_put/access_layer_language: Programming language of the Access Layer high level API used to PUT this IDS
* dataset_fair/ids_properties/provenance/node{:}/path: Path of the node within the IDS, following the syntax given in the link below. If empty, means the provenance information applies to the whole IDS.
* dataset_fair/ids_properties/provenance/node{:}/reference{:}/name: Reference name
* dataset_fair/ids_properties/provenance/node{:}/reference{:}/timestamp: Date and time (UTC) at which the reference was created, expressed in a human readable form (ISO 8601) : the format of the string shall be : YYYY-MM-DDTHH:MM:SSZ. Example : 2020-07-24T14:19:00Z


IDS summary
------------

* summary/time: Generic time
* summary/ids_properties/comment: Any comment describing the content of this IDS
* summary/ids_properties/homogeneous_time: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* summary/ids_properties/source: Source of the data (any comment describing the origin of the data : code, path to diagnostic signals, processing method, ...). Superseeded by the new provenance structure.
* summary/ids_properties/provider: Name of the person in charge of producing this data
* summary/ids_properties/creation_date: Date at which this data has been produced
* summary/code/name: Name of software generating IDS
* summary/code/commit: Unique commit reference of software
* summary/code/version: Unique version (tag) of software
* summary/code/repository: URL of software repository
* summary/code/parameters: List of the code specific parameters in XML format
* summary/code/output_flag: Output flag : 0 means the run is successful, other values mean some difficulty has been encountered, the exact meaning is then code specific. Negative values mean the result shall not be used.
* summary/global_quantities/ip/value: Value
* summary/global_quantities/current_non_inductive/value: Value
* summary/global_quantities/current_bootstrap/value: Value
* summary/global_quantities/current_ohm/value: Value
* summary/global_quantities/current_alignment/value: Value
* summary/global_quantities/v_loop/value: Value
* summary/global_quantities/li/value: Value
* summary/global_quantities/beta_tor/value: Value
* summary/global_quantities/beta_tor_norm/value: Value
* summary/global_quantities/beta_pol/value: Value
* summary/global_quantities/energy_diamagnetic/value: Value
* summary/global_quantities/energy_total/value: Value
* summary/global_quantities/energy_thermal/value: Value
* summary/global_quantities/energy_b_field_pol/value: Value
* summary/global_quantities/volume/value: Value
* summary/global_quantities/h_mode/value: Value
* summary/global_quantities/r0/value: Value
* summary/global_quantities/b0/value: Value
* summary/global_quantities/fusion_gain/value: Value
* summary/global_quantities/h_98/value: Value
* summary/global_quantities/tau_energy/value: Value
* summary/global_quantities/tau_helium/value: Value
* summary/global_quantities/tau_resistive/value: Value
* summary/global_quantities/resistance/value: Value
* summary/global_quantities/power_ohm/value: Value
* summary/global_quantities/power_steady/value: Value
* summary/global_quantities/power_radiated/value: Value
* summary/global_quantities/power_radiated_inside_lcfs/value: Value
* summary/global_quantities/power_radiated_outside_lcfs/value: Value
* summary/global_quantities/power_line/value: Value
* summary/global_quantities/power_bremsstrahlung/value: Value
* summary/global_quantities/power_synchrotron/value: Value
* summary/global_quantities/power_loss/value: Value
* summary/global_quantities/greenwald_fraction/value: Value
* summary/global_quantities/fusion_fluence/value: Value
* summary/global_quantities/psi_external_average/value: Value
* summary/line_average/t_e/value: Value
* summary/line_average/t_i_average/value: Value
* summary/line_average/n_e/value: Value
* summary/line_average/n_i_total/value: Value
* summary/line_average/zeff/value: Value
* summary/volume_average/t_e/value: Value
* summary/volume_average/n_e/value: Value
* summary/volume_average/n_i_total/value: Value
* summary/volume_average/zeff/value: Value
* summary/fusion/power/value: Value
* summary/fusion/current/value: Value
* summary/fusion/neutron_power_total/value: Value
* summary/runaways/particles/value: Value
* summary/runaways/current/value: Value
* summary/scrape_off_layer/t_e_decay_length/value: Value
* summary/scrape_off_layer/t_i_average_decay_length/value: Value
* summary/scrape_off_layer/n_e_decay_length/value: Value
* summary/scrape_off_layer/n_i_total_decay_length/value: Value
* summary/scrape_off_layer/heat_flux_e_decay_length/value: Value
* summary/scrape_off_layer/heat_flux_i_decay_length/value: Value
* summary/scrape_off_layer/power_radiated/value: Value
* summary/scrape_off_layer/pressure_neutral/value: Value
* summary/ids_properties/provenance/node{:}/path: Path of the node within the IDS, following the syntax given in the link below. If empty, means the provenance information applies to the whole IDS.
* summary/local/magnetic_axis/position/rho_tor_norm: Normalised toroidal flux coordinate. The normalizing value for rho_tor_norm, is the toroidal flux coordinate at the equilibrium boundary (LCFS or 99.x % of the LCFS in case of a fixed boundary equilibium calculation, see time_slice/boundary/b_flux_pol_norm in the equilibrium IDS) 
* summary/local/magnetic_axis/position/rho_tor: Toroidal flux coordinate. rho_tor = sqrt(b_flux_tor/(pi*b0)) ~ sqrt(pi*r^2*b0/(pi*b0)) ~ r [m]. The toroidal field used in its definition is indicated under global_quantities/b0
* summary/local/magnetic_axis/position/psi: Poloidal magnetic flux
* summary/local/magnetic_axis/position/r: Major radius
* summary/local/magnetic_axis/position/z: Height
* summary/local/magnetic_axis/t_e/value: Value
* summary/local/magnetic_axis/t_i_average/value: Value
* summary/local/magnetic_axis/n_e/value: Value
* summary/local/magnetic_axis/n_i_total/value: Value
* summary/local/magnetic_axis/zeff/value: Value
* summary/local/magnetic_axis/momentum_tor/value: Value
* summary/local/magnetic_axis/q/value: Value
* summary/local/magnetic_axis/magnetic_shear/value: Value
* summary/local/pedestal/position/rho_tor_norm: Normalised toroidal flux coordinate. The normalizing value for rho_tor_norm, is the toroidal flux coordinate at the equilibrium boundary (LCFS or 99.x % of the LCFS in case of a fixed boundary equilibium calculation, see time_slice/boundary/b_flux_pol_norm in the equilibrium IDS) 
* summary/local/pedestal/position/rho_tor: Toroidal flux coordinate. rho_tor = sqrt(b_flux_tor/(pi*b0)) ~ sqrt(pi*r^2*b0/(pi*b0)) ~ r [m]. The toroidal field used in its definition is indicated under global_quantities/b0
* summary/local/pedestal/position/psi: Poloidal magnetic flux
* summary/local/pedestal/t_i_average/value: Value
* summary/local/pedestal/n_e/value: Value
* summary/local/pedestal/n_i_total/value: Value
* summary/local/pedestal/zeff/value: Value
* summary/local/pedestal/momentum_tor/value: Value
* summary/local/pedestal/q/value: Value
* summary/local/pedestal/magnetic_shear/value: Value
* summary/local/itb/position/rho_tor_norm: Normalised toroidal flux coordinate. The normalizing value for rho_tor_norm, is the toroidal flux coordinate at the equilibrium boundary (LCFS or 99.x % of the LCFS in case of a fixed boundary equilibium calculation, see time_slice/boundary/b_flux_pol_norm in the equilibrium IDS) 
* summary/local/itb/position/rho_tor: Toroidal flux coordinate. rho_tor = sqrt(b_flux_tor/(pi*b0)) ~ sqrt(pi*r^2*b0/(pi*b0)) ~ r [m]. The toroidal field used in its definition is indicated under global_quantities/b0
* summary/local/itb/position/psi: Poloidal magnetic flux
* summary/local/itb/t_e/value: Value
* summary/local/itb/t_i_average/value: Value
* summary/local/itb/n_e/value: Value
* summary/local/itb/n_i_total/value: Value
* summary/local/itb/zeff/value: Value
* summary/local/itb/momentum_tor/value: Value
* summary/local/itb/q/value: Value
* summary/local/itb/magnetic_shear/value: Value
* summary/local/limiter/t_e/value: Value
* summary/local/limiter/n_e/value: Value
* summary/local/divertor_plate{:}/t_e/value: Value
* summary/local/divertor_plate{:}/n_e/value: Value
* summary/local/divertor_plate{:}/zeff/value: Value
* summary/line_average/n_i/hydrogen/value: Value
* summary/line_average/n_i/deuterium/value: Value
* summary/line_average/n_i/tritium/value: Value
* summary/line_average/n_i/helium_3/value: Value
* summary/line_average/n_i/helium_4/value: Value
* summary/line_average/n_i/lithium/value: Value
* summary/line_average/n_i/carbon/value: Value
* summary/line_average/n_i/nitrogen/value: Value
* summary/line_average/n_i/neon/value: Value
* summary/line_average/n_i/argon/value: Value
* summary/line_average/n_i/xenon/value: Value
* summary/line_average/n_i/oxygen/value: Value
* summary/line_average/n_i/tungsten/value: Value
* summary/volume_average/n_i/hydrogen/value: Value
* summary/volume_average/n_i/deuterium/value: Value
* summary/volume_average/n_i/tritium/value: Value
* summary/volume_average/n_i/helium_3/value: Value
* summary/volume_average/n_i/helium_4/value: Value
* summary/volume_average/n_i/lithium/value: Value
* summary/volume_average/n_i/carbon/value: Value
* summary/volume_average/n_i/nitrogen/value: Value
* summary/volume_average/n_i/neon/value: Value
* summary/volume_average/n_i/argon/value: Value
* summary/volume_average/n_i/xenon/value: Value
* summary/volume_average/n_i/oxygen/value: Value
* summary/volume_average/n_i/tungsten/value: Value
* summary/fusion/neutron_fluxes/total/value: Value
* summary/fusion/neutron_fluxes/thermal/value: Value
* summary/fusion/neutron_rates/total/value: Value
* summary/fusion/neutron_rates/thermal/value: Value
* summary/heating_current_drive/ec{:}/position/value: Value
* summary/heating_current_drive/ec{:}/power/value: Value
* summary/heating_current_drive/ec{:}/current/value: Value
* summary/heating_current_drive/nbi{:}/power/value: Value
* summary/heating_current_drive/nbi{:}/current/value: Value
* summary/heating_current_drive/nbi{:}/tangency_radius/value: Value
* summary/heating_current_drive/nbi{:}/angle/value: Value
* summary/heating_current_drive/nbi{:}/direction/value: Value
* summary/heating_current_drive/nbi{:}/energy/value: Value
* summary/heating_current_drive/ic{:}/frequency/value: Value
* summary/heating_current_drive/ic{:}/position/value: Value
* summary/heating_current_drive/ic{:}/n_tor/value: Value
* summary/heating_current_drive/ic{:}/harmonic/value: Value
* summary/heating_current_drive/ic{:}/power/value: Value
* summary/heating_current_drive/ic{:}/current/value: Value
* summary/heating_current_drive/ic{:}/energy_fast/value: Value
* summary/heating_current_drive/lh{:}/frequency/value: Value
* summary/heating_current_drive/lh{:}/position/value: Value
* summary/heating_current_drive/lh{:}/n_parallel/value: Value
* summary/heating_current_drive/lh{:}/power/value: Value
* summary/heating_current_drive/lh{:}/current/value: Value
* summary/heating_current_drive/lh{:}/energy_fast/value: Value
* summary/ids_properties/provenance/node{:}/reference{:}/name: Reference name
* summary/ids_properties/provenance/node{:}/reference{:}/timestamp: Date and time (UTC) at which the reference was created, expressed in a human readable form (ISO 8601) : the format of the string shall be : YYYY-MM-DDTHH:MM:SSZ. Example : 2020-07-24T14:19:00Z
* summary/local/magnetic_axis/n_i/hydrogen/value: Value
* summary/local/magnetic_axis/n_i/deuterium/value: Value
* summary/local/magnetic_axis/n_i/tritium/value: Value
* summary/local/magnetic_axis/n_i/helium_3/value: Value
* summary/local/magnetic_axis/n_i/helium_4/value: Value
* summary/local/magnetic_axis/n_i/lithium/value: Value
* summary/local/magnetic_axis/n_i/carbon/value: Value
* summary/local/magnetic_axis/n_i/nitrogen/value: Value
* summary/local/magnetic_axis/n_i/neon/value: Value
* summary/local/magnetic_axis/n_i/argon/value: Value
* summary/local/magnetic_axis/n_i/xenon/value: Value
* summary/local/magnetic_axis/n_i/oxygen/value: Value
* summary/local/magnetic_axis/n_i/tungsten/value: Value
* summary/local/magnetic_axis/velocity_tor/hydrogen/value: Value
* summary/local/magnetic_axis/velocity_tor/deuterium/value: Value
* summary/local/magnetic_axis/velocity_tor/tritium/value: Value
* summary/local/magnetic_axis/velocity_tor/helium_3/value: Value
* summary/local/magnetic_axis/velocity_tor/helium_4/value: Value
* summary/local/magnetic_axis/velocity_tor/lithium/value: Value
* summary/local/magnetic_axis/velocity_tor/carbon/value: Value
* summary/local/magnetic_axis/velocity_tor/nitrogen/value: Value
* summary/local/magnetic_axis/velocity_tor/neon/value: Value
* summary/local/magnetic_axis/velocity_tor/argon/value: Value
* summary/local/magnetic_axis/velocity_tor/xenon/value: Value
* summary/local/magnetic_axis/velocity_tor/oxygen/value: Value
* summary/local/magnetic_axis/velocity_tor/tungsten/value: Value
* summary/local/magnetic_axis/velocity_phi/hydrogen/value: Value
* summary/local/magnetic_axis/velocity_phi/deuterium/value: Value
* summary/local/magnetic_axis/velocity_phi/tritium/value: Value
* summary/local/magnetic_axis/velocity_phi/helium_3/value: Value
* summary/local/magnetic_axis/velocity_phi/helium_4/value: Value
* summary/local/magnetic_axis/velocity_phi/lithium/value: Value
* summary/local/magnetic_axis/velocity_phi/carbon/value: Value
* summary/local/magnetic_axis/velocity_phi/nitrogen/value: Value
* summary/local/magnetic_axis/velocity_phi/neon/value: Value
* summary/local/magnetic_axis/velocity_phi/argon/value: Value
* summary/local/magnetic_axis/velocity_phi/xenon/value: Value
* summary/local/magnetic_axis/velocity_phi/oxygen/value: Value
* summary/local/magnetic_axis/velocity_phi/tungsten/value: Value
* summary/local/pedestal/n_i/hydrogen/value: Value
* summary/local/pedestal/n_i/deuterium/value: Value
* summary/local/pedestal/n_i/tritium/value: Value
* summary/local/pedestal/n_i/helium_3/value: Value
* summary/local/pedestal/n_i/helium_4/value: Value
* summary/local/pedestal/n_i/lithium/value: Value
* summary/local/pedestal/n_i/carbon/value: Value
* summary/local/pedestal/n_i/nitrogen/value: Value
* summary/local/pedestal/n_i/neon/value: Value
* summary/local/pedestal/n_i/argon/value: Value
* summary/local/pedestal/n_i/xenon/value: Value
* summary/local/pedestal/n_i/oxygen/value: Value
* summary/local/pedestal/n_i/tungsten/value: Value
* summary/local/pedestal/velocity_tor/hydrogen/value: Value
* summary/local/pedestal/velocity_tor/deuterium/value: Value
* summary/local/pedestal/velocity_tor/tritium/value: Value
* summary/local/pedestal/velocity_tor/helium_3/value: Value
* summary/local/pedestal/velocity_tor/helium_4/value: Value
* summary/local/pedestal/velocity_tor/lithium/value: Value
* summary/local/pedestal/velocity_tor/carbon/value: Value
* summary/local/pedestal/velocity_tor/nitrogen/value: Value
* summary/local/pedestal/velocity_tor/neon/value: Value
* summary/local/pedestal/velocity_tor/argon/value: Value
* summary/local/pedestal/velocity_tor/xenon/value: Value
* summary/local/pedestal/velocity_tor/oxygen/value: Value
* summary/local/pedestal/velocity_tor/tungsten/value: Value
* summary/local/pedestal/velocity_phi/hydrogen/value: Value
* summary/local/pedestal/velocity_phi/deuterium/value: Value
* summary/local/pedestal/velocity_phi/tritium/value: Value
* summary/local/pedestal/velocity_phi/helium_3/value: Value
* summary/local/pedestal/velocity_phi/helium_4/value: Value
* summary/local/pedestal/velocity_phi/lithium/value: Value
* summary/local/pedestal/velocity_phi/carbon/value: Value
* summary/local/pedestal/velocity_phi/nitrogen/value: Value
* summary/local/pedestal/velocity_phi/neon/value: Value
* summary/local/pedestal/velocity_phi/argon/value: Value
* summary/local/pedestal/velocity_phi/xenon/value: Value
* summary/local/pedestal/velocity_phi/oxygen/value: Value
* summary/local/pedestal/velocity_phi/tungsten/value: Value
* summary/local/itb/n_i/hydrogen/value: Value
* summary/local/itb/n_i/deuterium/value: Value
* summary/local/itb/n_i/tritium/value: Value
* summary/local/itb/n_i/helium_3/value: Value
* summary/local/itb/n_i/helium_4/value: Value
* summary/local/itb/n_i/lithium/value: Value
* summary/local/itb/n_i/carbon/value: Value
* summary/local/itb/n_i/nitrogen/value: Value
* summary/local/itb/n_i/neon/value: Value
* summary/local/itb/n_i/argon/value: Value
* summary/local/itb/n_i/xenon/value: Value
* summary/local/itb/n_i/oxygen/value: Value
* summary/local/itb/n_i/tungsten/value: Value
* summary/local/itb/velocity_tor/hydrogen/value: Value
* summary/local/itb/velocity_tor/deuterium/value: Value
* summary/local/itb/velocity_tor/tritium/value: Value
* summary/local/itb/velocity_tor/helium_3/value: Value
* summary/local/itb/velocity_tor/helium_4/value: Value
* summary/local/itb/velocity_tor/lithium/value: Value
* summary/local/itb/velocity_tor/carbon/value: Value
* summary/local/itb/velocity_tor/nitrogen/value: Value
* summary/local/itb/velocity_tor/neon/value: Value
* summary/local/itb/velocity_tor/argon/value: Value
* summary/local/itb/velocity_tor/xenon/value: Value
* summary/local/itb/velocity_tor/oxygen/value: Value
* summary/local/itb/velocity_tor/tungsten/value: Value
* summary/local/itb/velocity_phi/hydrogen/value: Value
* summary/local/itb/velocity_phi/deuterium/value: Value
* summary/local/itb/velocity_phi/tritium/value: Value
* summary/local/itb/velocity_phi/helium_3/value: Value
* summary/local/itb/velocity_phi/helium_4/value: Value
* summary/local/itb/velocity_phi/lithium/value: Value
* summary/local/itb/velocity_phi/carbon/value: Value
* summary/local/itb/velocity_phi/nitrogen/value: Value
* summary/local/itb/velocity_phi/neon/value: Value
* summary/local/itb/velocity_phi/argon/value: Value
* summary/local/itb/velocity_phi/xenon/value: Value
* summary/local/itb/velocity_phi/oxygen/value: Value
* summary/local/itb/velocity_phi/tungsten/value: Value
* summary/fusion/neutron_fluxes/dd/total/value: Value
* summary/fusion/neutron_fluxes/dd/thermal/value: Value
* summary/fusion/neutron_fluxes/dd/beam_thermal/value: Value
* summary/fusion/neutron_fluxes/dd/beam_beam/value: Value
* summary/fusion/neutron_fluxes/dt/total/value: Value
* summary/fusion/neutron_fluxes/dt/thermal/value: Value
* summary/fusion/neutron_fluxes/dt/beam_thermal/value: Value
* summary/fusion/neutron_fluxes/dt/beam_beam/value: Value
* summary/fusion/neutron_fluxes/tt/total/value: Value
* summary/fusion/neutron_fluxes/tt/thermal/value: Value
* summary/fusion/neutron_fluxes/tt/beam_thermal/value: Value
* summary/fusion/neutron_fluxes/tt/beam_beam/value: Value
* summary/fusion/neutron_rates/dd/total/value: Value
* summary/fusion/neutron_rates/dd/thermal/value: Value
* summary/fusion/neutron_rates/dd/beam_thermal/value: Value
* summary/fusion/neutron_rates/dd/beam_beam/value: Value
* summary/fusion/neutron_rates/dt/total/value: Value
* summary/fusion/neutron_rates/dt/thermal/value: Value
* summary/fusion/neutron_rates/dt/beam_thermal/value: Value
* summary/fusion/neutron_rates/dt/beam_beam/value: Value
* summary/fusion/neutron_rates/tt/total/value: Value
* summary/fusion/neutron_rates/tt/thermal/value: Value
* summary/fusion/neutron_rates/tt/beam_thermal/value: Value
* summary/fusion/neutron_rates/tt/beam_beam/value: Value
* summary/heating_current_drive/nbi{:}/species/a/value: Value
* summary/heating_current_drive/nbi{:}/species/z_n/value: Value


IDS core_profiles
------------------

* core_profiles/time: Generic time
* core_profiles/ids_properties/comment: Any comment describing the content of this IDS
* core_profiles/ids_properties/homogeneous_time: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* core_profiles/ids_properties/source: Source of the data (any comment describing the origin of the data : code, path to diagnostic signals, processing method, ...). Superseeded by the new provenance structure.
* core_profiles/ids_properties/provider: Name of the person in charge of producing this data
* core_profiles/ids_properties/creation_date: Date at which this data has been produced
* core_profiles/profiles_1d{:}/t_i_average: Ion temperature (averaged on charge states and ion species)
* core_profiles/profiles_1d{:}/n_i_total_over_n_e: Ratio of total ion density (sum over species and charge states) over electron density. (thermal+non-thermal)
* core_profiles/profiles_1d{:}/momentum_tor: Total plasma toroidal momentum, summed over ion species and electrons weighted by their density and major radius, i.e. sum_over_species(n*R*m*Vphi)
* core_profiles/profiles_1d{:}/momentum_phi: Total plasma toroidal momentum, summed over ion species and electrons weighted by their density and major radius, i.e. sum_over_species(n*R*m*Vphi)
* core_profiles/profiles_1d{:}/zeff: Effective charge
* core_profiles/profiles_1d{:}/pressure_ion_total: Total (sum over ion species) thermal ion pressure
* core_profiles/profiles_1d{:}/pressure_thermal: Thermal pressure (electrons+ions)
* core_profiles/profiles_1d{:}/pressure_perpendicular: Total perpendicular pressure (electrons+ions, thermal+non-thermal)
* core_profiles/profiles_1d{:}/pressure_parallel: Total parallel pressure (electrons+ions, thermal+non-thermal)
* core_profiles/profiles_1d{:}/j_total: Total parallel current density = average(jtot.B) / B0, where B0 = Core_Profiles/Vacuum_Toroidal_Field/ B0
* core_profiles/profiles_1d{:}/current_parallel_inside: Parallel current driven inside the flux surface. Cumulative surface integral of j_total
* core_profiles/profiles_1d{:}/j_tor: Total toroidal current density = average(J_Tor/R) / average(1/R)
* core_profiles/profiles_1d{:}/j_phi: Total toroidal current density = average(J_phi/R) / average(1/R)
* core_profiles/profiles_1d{:}/j_ohmic: Ohmic parallel current density = average(J_Ohmic.B) / B0, where B0 = Core_Profiles/Vacuum_Toroidal_Field/ B0
* core_profiles/profiles_1d{:}/j_non_inductive: Non-inductive (includes bootstrap) parallel current density = average(jni.B) / B0, where B0 = Core_Profiles/Vacuum_Toroidal_Field/ B0
* core_profiles/profiles_1d{:}/j_bootstrap: Bootstrap current density = average(J_Bootstrap.B) / B0, where B0 = Core_Profiles/Vacuum_Toroidal_Field/ B0
* core_profiles/profiles_1d{:}/conductivity_parallel: Parallel conductivity
* core_profiles/profiles_1d{:}/e_field_parallel: Parallel electric field = average(E.B) / B0, where Core_Profiles/Vacuum_Toroidal_Field/ B0
* core_profiles/profiles_1d{:}/q: Safety factor (IMAS uses COCOS=11: only positive when toroidal current and magnetic field are in same direction)
* core_profiles/profiles_1d{:}/magnetic_shear: Magnetic shear, defined as rho_tor/q . dq/drho_tor
* core_profiles/profiles_1d{:}/time: Time
* core_profiles/global_quantities/ip: Total plasma current (toroidal component). Positive sign means anti-clockwise when viewed from above.
* core_profiles/global_quantities/current_non_inductive: Total non-inductive current (toroidal component). Positive sign means anti-clockwise when viewed from above.
* core_profiles/global_quantities/current_bootstrap: Bootstrap current (toroidal component). Positive sign means anti-clockwise when viewed from above.
* core_profiles/global_quantities/v_loop: LCFS loop voltage (positive value drives positive ohmic current that flows anti-clockwise when viewed from above)
* core_profiles/global_quantities/li: Internal inductance. The li_3 definition is used, i.e. li_3 = 2/R0/mu0^2/Ip^2 * int(Bp^2 dV).
* core_profiles/global_quantities/li_3: Internal inductance. The li_3 definition is used, i.e. li_3 = 2/R0/mu0^2/Ip^2 * int(Bp^2 dV).
* core_profiles/global_quantities/beta_tor: Toroidal beta, defined as the volume-averaged total perpendicular pressure divided by (B0^2/(2*mu0)), i.e. beta_toroidal = 2 mu0 int(p dV) / V / B0^2
* core_profiles/global_quantities/beta_tor_norm: Normalised toroidal beta, defined as 100 * beta_tor * a[m] * B0 [T] / ip [MA] 
* core_profiles/global_quantities/beta_pol: Poloidal beta. Defined as betap = 4 int(p dV) / [R_0 * mu_0 * Ip^2]
* core_profiles/global_quantities/energy_diamagnetic: Plasma energy content = 3/2 * integral over the plasma volume of the total perpendicular pressure 
* core_profiles/vacuum_toroidal_field/r0: Reference major radius where the vacuum toroidal magnetic field is given (usually a fixed position such as the middle of the vessel at the equatorial midplane)
* core_profiles/vacuum_toroidal_field/b0: Vacuum toroidal field at R0 [T]; Positive sign means anti-clockwise when viewing from above. The product R0B0 must be consistent with the b_tor_vacuum_r field of the tf IDS.
* core_profiles/code/name: Name of software generating IDS
* core_profiles/code/commit: Unique commit reference of software
* core_profiles/code/version: Unique version (tag) of software
* core_profiles/code/repository: URL of software repository
* core_profiles/code/parameters: List of the code specific parameters in XML format
* core_profiles/code/output_flag: Output flag : 0 means the run is successful, other values mean some difficulty has been encountered, the exact meaning is then code specific. Negative values mean the result shall not be used.
* core_profiles/profiles_1d{:}/grid/rho_tor_norm: Normalised toroidal flux coordinate. The normalizing value for rho_tor_norm, is the toroidal flux coordinate at the equilibrium boundary (LCFS or 99.x % of the LCFS in case of a fixed boundary equilibium calculation, see time_slice/boundary/b_flux_pol_norm in the equilibrium IDS) 
* core_profiles/profiles_1d{:}/grid/rho_tor: Toroidal flux coordinate = sqrt(phi/(pi*b0)), where the toroidal magnetic field, b0, corresponds to that stored in vacuum_toroidal_field/b0 and pi can be found in the IMAS constants
* core_profiles/profiles_1d{:}/grid/rho_pol_norm: Normalised poloidal flux coordinate = sqrt((psi(rho)-psi(magnetic_axis)) / (psi(LCFS)-psi(magnetic_axis)))
* core_profiles/profiles_1d{:}/grid/psi: Poloidal magnetic flux
* core_profiles/profiles_1d{:}/grid/volume: Volume enclosed inside the magnetic surface
* core_profiles/profiles_1d{:}/grid/area: Cross-sectional area of the flux surface
* core_profiles/profiles_1d{:}/grid/surface: Surface area of the toroidal flux surface
* core_profiles/profiles_1d{:}/grid/psi_magnetic_axis: Value of the poloidal magnetic flux at the magnetic axis (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* core_profiles/profiles_1d{:}/grid/psi_boundary: Value of the poloidal magnetic flux at the plasma boundary (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* core_profiles/profiles_1d{:}/electrons/temperature: Temperature
* core_profiles/profiles_1d{:}/electrons/density: Density (thermal+non-thermal)
* core_profiles/profiles_1d{:}/electrons/density_thermal: Density of thermal particles
* core_profiles/profiles_1d{:}/electrons/density_fast: Density of fast (non-thermal) particles
* core_profiles/profiles_1d{:}/electrons/pressure: Pressure (thermal+non-thermal)
* core_profiles/profiles_1d{:}/electrons/pressure_fast_perpendicular: Fast (non-thermal) perpendicular pressure
* core_profiles/profiles_1d{:}/electrons/pressure_fast_parallel: Fast (non-thermal) parallel pressure
* core_profiles/profiles_1d{:}/ion{:}/z_ion: Ion charge (of the dominant ionisation state; lumped ions are allowed), volume averaged over plasma radius
* core_profiles/profiles_1d{:}/ion{:}/label: String identifying ion (e.g. H, D, T, He, C, D2, ...)
* core_profiles/profiles_1d{:}/ion{:}/name: String identifying ion (e.g. H, D, T, He, C, D2, ...)
* core_profiles/profiles_1d{:}/ion{:}/temperature: Temperature (average over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/ion{:}/density: Density (thermal+non-thermal) (sum over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/ion{:}/density_thermal: Density (thermal) (sum over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/ion{:}/density_fast: Density of fast (non-thermal) particles (sum over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/ion{:}/pressure: Pressure (thermal+non-thermal) (sum over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/ion{:}/pressure_fast_perpendicular: Fast (non-thermal) perpendicular pressure  (sum over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/ion{:}/pressure_fast_parallel: Fast (non-thermal) parallel pressure  (sum over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/ion{:}/velocity_tor: Toroidal velocity (average over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/ion{:}/velocity_pol: Poloidal velocity (average over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/neutral{:}/label: String identifying the species (e.g. H, D, T, He, C, D2, DT, CD4, ...)
* core_profiles/profiles_1d{:}/neutral{:}/name: String identifying the species (e.g. H, D, T, He, C, D2, DT, CD4, ...)
* core_profiles/profiles_1d{:}/neutral{:}/ion_index: Index of the corresponding ion species in the ../../ion array
* core_profiles/profiles_1d{:}/neutral{:}/temperature: Temperature (average over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/neutral{:}/density: Density (thermal+non-thermal) (sum over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/neutral{:}/density_thermal: Density (thermal) (sum over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/neutral{:}/density_fast: Density of fast (non-thermal) particles (sum over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/neutral{:}/pressure: Pressure (thermal+non-thermal) (sum over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/neutral{:}/pressure_fast_perpendicular: Fast (non-thermal) perpendicular pressure  (sum over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/neutral{:}/pressure_fast_parallel: Fast (non-thermal) parallel pressure  (sum over charge states when multiple charge states are considered)
* core_profiles/profiles_1d{:}/e_field/radial: Radial component
* core_profiles/profiles_1d{:}/e_field/parallel: Parallel component
* core_profiles/ids_properties/provenance/node{:}/path: Path of the node within the IDS, following the syntax given in the link below. If empty, means the provenance information applies to the whole IDS.
* core_profiles/profiles_1d{:}/ion{:}/element{:}/a: Mass of atom
* core_profiles/profiles_1d{:}/ion{:}/element{:}/z_n: Nuclear charge
* core_profiles/profiles_1d{:}/ion{:}/temperature_fit/source: Path to the source data for each measurement in the IMAS data dictionary, e.g. ece/channel(i)/t_e for the electron temperature on the i-th channel in the ECE IDS
* core_profiles/profiles_1d{:}/ion{:}/density_fit/source: Path to the source data for each measurement in the IMAS data dictionary, e.g. ece/channel(i)/t_e for the electron temperature on the i-th channel in the ECE IDS
* core_profiles/profiles_1d{:}/ion{:}/velocity/poloidal: Poloidal component
* core_profiles/profiles_1d{:}/ion{:}/velocity/toroidal: Toroidal component
* core_profiles/profiles_1d{:}/neutral{:}/element{:}/a: Mass of atom
* core_profiles/profiles_1d{:}/neutral{:}/element{:}/z_n: Nuclear charge
* core_profiles/ids_properties/provenance/node{:}/reference{:}/name: Reference name
* core_profiles/ids_properties/provenance/node{:}/reference{:}/timestamp: Date and time (UTC) at which the reference was created, expressed in a human readable form (ISO 8601) : the format of the string shall be : YYYY-MM-DDTHH:MM:SSZ. Example : 2020-07-24T14:19:00Z
* core_profiles/profiles_1d{:}/ion{:}/state{:}/density_fit/source: Path to the source data for each measurement in the IMAS data dictionary, e.g. ece/channel(i)/t_e for the electron temperature on the i-th channel in the ECE IDS


IDS core_transport
-------------------

* core_transport/time: Generic time
* core_transport/ids_properties/comment: Any comment describing the content of this IDS
* core_transport/ids_properties/homogeneous_time: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* core_transport/ids_properties/source: Source of the data (any comment describing the origin of the data : code, path to diagnostic signals, processing method, ...). Superseeded by the new provenance structure.
* core_transport/ids_properties/provider: Name of the person in charge of producing this data
* core_transport/ids_properties/creation_date: Date at which this data has been produced
* core_transport/vacuum_toroidal_field/r0: Reference major radius where the vacuum toroidal magnetic field is given (usually a fixed position such as the middle of the vessel at the equatorial midplane)
* core_transport/vacuum_toroidal_field/b0: Vacuum toroidal field at R0 [T]; Positive sign means anti-clockwise when viewing from above. The product R0B0 must be consistent with the b_tor_vacuum_r field of the tf IDS.
* core_transport/model{:}/flux_multiplier: Multiplier applied to the particule flux when adding its contribution in the expression of the heat flux : can be 0, 3/2 or 5/2
* core_transport/code/name: Name of software generating IDS
* core_transport/code/commit: Unique commit reference of software
* core_transport/code/version: Unique version (tag) of software
* core_transport/code/repository: URL of software repository
* core_transport/code/parameters: List of the code specific parameters in XML format
* core_transport/code/output_flag: Output flag : 0 means the run is successful, other values mean some difficulty has been encountered, the exact meaning is then code specific. Negative values mean the result shall not be used.
* core_transport/model{:}/identifier/name: Short string identifier
* core_transport/model{:}/identifier/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* core_transport/model{:}/identifier/description: Verbose description
* core_transport/model{:}/profiles_1d{:}/conductivity_parallel: Parallel conductivity
* core_transport/model{:}/profiles_1d{:}/time: Time
* core_transport/ids_properties/provenance/node{:}/path: Path of the node within the IDS, following the syntax given in the link below. If empty, means the provenance information applies to the whole IDS.
* core_transport/model{:}/profiles_1d{:}/grid_d/rho_tor_norm: Normalised toroidal flux coordinate. The normalizing value for rho_tor_norm, is the toroidal flux coordinate at the equilibrium boundary (LCFS or 99.x % of the LCFS in case of a fixed boundary equilibium calculation, see time_slice/boundary/b_flux_pol_norm in the equilibrium IDS) 
* core_transport/model{:}/profiles_1d{:}/grid_d/rho_tor: Toroidal flux coordinate = sqrt(phi/(pi*b0)), where the toroidal magnetic field, b0, corresponds to that stored in vacuum_toroidal_field/b0 and pi can be found in the IMAS constants
* core_transport/model{:}/profiles_1d{:}/grid_d/rho_pol_norm: Normalised poloidal flux coordinate = sqrt((psi(rho)-psi(magnetic_axis)) / (psi(LCFS)-psi(magnetic_axis)))
* core_transport/model{:}/profiles_1d{:}/grid_d/psi: Poloidal magnetic flux
* core_transport/model{:}/profiles_1d{:}/grid_d/volume: Volume enclosed inside the magnetic surface
* core_transport/model{:}/profiles_1d{:}/grid_d/area: Cross-sectional area of the flux surface
* core_transport/model{:}/profiles_1d{:}/grid_d/surface: Surface area of the toroidal flux surface
* core_transport/model{:}/profiles_1d{:}/grid_d/psi_magnetic_axis: Value of the poloidal magnetic flux at the magnetic axis (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* core_transport/model{:}/profiles_1d{:}/grid_d/psi_boundary: Value of the poloidal magnetic flux at the plasma boundary (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* core_transport/model{:}/profiles_1d{:}/grid_v/rho_tor_norm: Normalised toroidal flux coordinate. The normalizing value for rho_tor_norm, is the toroidal flux coordinate at the equilibrium boundary (LCFS or 99.x % of the LCFS in case of a fixed boundary equilibium calculation, see time_slice/boundary/b_flux_pol_norm in the equilibrium IDS) 
* core_transport/model{:}/profiles_1d{:}/grid_v/rho_tor: Toroidal flux coordinate = sqrt(phi/(pi*b0)), where the toroidal magnetic field, b0, corresponds to that stored in vacuum_toroidal_field/b0 and pi can be found in the IMAS constants
* core_transport/model{:}/profiles_1d{:}/grid_v/rho_pol_norm: Normalised poloidal flux coordinate = sqrt((psi(rho)-psi(magnetic_axis)) / (psi(LCFS)-psi(magnetic_axis)))
* core_transport/model{:}/profiles_1d{:}/grid_v/psi: Poloidal magnetic flux
* core_transport/model{:}/profiles_1d{:}/grid_v/volume: Volume enclosed inside the magnetic surface
* core_transport/model{:}/profiles_1d{:}/grid_v/area: Cross-sectional area of the flux surface
* core_transport/model{:}/profiles_1d{:}/grid_v/surface: Surface area of the toroidal flux surface
* core_transport/model{:}/profiles_1d{:}/grid_v/psi_magnetic_axis: Value of the poloidal magnetic flux at the magnetic axis (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* core_transport/model{:}/profiles_1d{:}/grid_v/psi_boundary: Value of the poloidal magnetic flux at the plasma boundary (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* core_transport/model{:}/profiles_1d{:}/grid_flux/rho_tor_norm: Normalised toroidal flux coordinate. The normalizing value for rho_tor_norm, is the toroidal flux coordinate at the equilibrium boundary (LCFS or 99.x % of the LCFS in case of a fixed boundary equilibium calculation, see time_slice/boundary/b_flux_pol_norm in the equilibrium IDS) 
* core_transport/model{:}/profiles_1d{:}/grid_flux/rho_tor: Toroidal flux coordinate = sqrt(phi/(pi*b0)), where the toroidal magnetic field, b0, corresponds to that stored in vacuum_toroidal_field/b0 and pi can be found in the IMAS constants
* core_transport/model{:}/profiles_1d{:}/grid_flux/rho_pol_norm: Normalised poloidal flux coordinate = sqrt((psi(rho)-psi(magnetic_axis)) / (psi(LCFS)-psi(magnetic_axis)))
* core_transport/model{:}/profiles_1d{:}/grid_flux/psi: Poloidal magnetic flux
* core_transport/model{:}/profiles_1d{:}/grid_flux/volume: Volume enclosed inside the magnetic surface
* core_transport/model{:}/profiles_1d{:}/grid_flux/area: Cross-sectional area of the flux surface
* core_transport/model{:}/profiles_1d{:}/grid_flux/surface: Surface area of the toroidal flux surface
* core_transport/model{:}/profiles_1d{:}/grid_flux/psi_magnetic_axis: Value of the poloidal magnetic flux at the magnetic axis (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* core_transport/model{:}/profiles_1d{:}/grid_flux/psi_boundary: Value of the poloidal magnetic flux at the plasma boundary (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* core_transport/model{:}/profiles_1d{:}/total_ion_energy/d: Effective diffusivity
* core_transport/model{:}/profiles_1d{:}/total_ion_energy/v: Effective convection
* core_transport/model{:}/profiles_1d{:}/total_ion_energy/flux: Flux
* core_transport/model{:}/profiles_1d{:}/momentum_tor/d: Effective diffusivity
* core_transport/model{:}/profiles_1d{:}/momentum_tor/v: Effective convection
* core_transport/model{:}/profiles_1d{:}/momentum_tor/flux: Flux
* core_transport/ids_properties/provenance/node{:}/reference{:}/name: Reference name
* core_transport/ids_properties/provenance/node{:}/reference{:}/timestamp: Date and time (UTC) at which the reference was created, expressed in a human readable form (ISO 8601) : the format of the string shall be : YYYY-MM-DDTHH:MM:SSZ. Example : 2020-07-24T14:19:00Z
* core_transport/model{:}/profiles_1d{:}/electrons/particles/d: Effective diffusivity
* core_transport/model{:}/profiles_1d{:}/electrons/particles/v: Effective convection
* core_transport/model{:}/profiles_1d{:}/electrons/particles/flux: Flux
* core_transport/model{:}/profiles_1d{:}/electrons/energy/d: Effective diffusivity
* core_transport/model{:}/profiles_1d{:}/electrons/energy/v: Effective convection
* core_transport/model{:}/profiles_1d{:}/electrons/energy/flux: Flux


IDS core_sources
-----------------

* core_sources/time: Generic time
* core_sources/ids_properties/comment: Any comment describing the content of this IDS
* core_sources/ids_properties/homogeneous_time: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* core_sources/ids_properties/source: Source of the data (any comment describing the origin of the data : code, path to diagnostic signals, processing method, ...). Superseeded by the new provenance structure.
* core_sources/ids_properties/provider: Name of the person in charge of producing this data
* core_sources/ids_properties/creation_date: Date at which this data has been produced
* core_sources/vacuum_toroidal_field/r0: Reference major radius where the vacuum toroidal magnetic field is given (usually a fixed position such as the middle of the vessel at the equatorial midplane)
* core_sources/vacuum_toroidal_field/b0: Vacuum toroidal field at R0 [T]; Positive sign means anti-clockwise when viewing from above. The product R0B0 must be consistent with the b_tor_vacuum_r field of the tf IDS.
* core_sources/code/name: Name of software generating IDS
* core_sources/code/commit: Unique commit reference of software
* core_sources/code/version: Unique version (tag) of software
* core_sources/code/repository: URL of software repository
* core_sources/code/parameters: List of the code specific parameters in XML format
* core_sources/code/output_flag: Output flag : 0 means the run is successful, other values mean some difficulty has been encountered, the exact meaning is then code specific. Negative values mean the result shall not be used.
* core_sources/source{:}/identifier/name: Short string identifier
* core_sources/source{:}/identifier/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* core_sources/source{:}/identifier/description: Verbose description
* core_sources/source{:}/global_quantities{:}/power: Total power coupled to the plasma
* core_sources/source{:}/global_quantities{:}/total_ion_power: Total power coupled to ion species (summed over ion  species)
* core_sources/source{:}/global_quantities{:}/torque_tor: Toroidal torque
* core_sources/source{:}/global_quantities{:}/current_parallel: Parallel current driven
* core_sources/source{:}/global_quantities{:}/time: Time
* core_sources/source{:}/profiles_1d{:}/total_ion_energy: Source term for the total (summed over ion  species) energy equation
* core_sources/source{:}/profiles_1d{:}/total_ion_power_inside: Total power coupled to ion species (summed over ion  species) inside the flux surface. Cumulative volume integral of the source term for the total ion energy equation
* core_sources/source{:}/profiles_1d{:}/momentum_tor: Source term for total toroidal momentum equation
* core_sources/source{:}/profiles_1d{:}/torque_tor_inside: Toroidal torque inside the flux surface. Cumulative volume integral of the source term for the total toroidal momentum equation
* core_sources/source{:}/profiles_1d{:}/j_parallel: Parallel current density source, average(J.B) / B0, where B0 = core_sources/vacuum_toroidal_field/b0 
* core_sources/source{:}/profiles_1d{:}/current_parallel_inside: Parallel current driven inside the flux surface. Cumulative surface integral of j_parallel
* core_sources/source{:}/profiles_1d{:}/conductivity_parallel: Parallel conductivity due to this source
* core_sources/source{:}/profiles_1d{:}/time: Time
* core_sources/ids_properties/provenance/node{:}/path: Path of the node within the IDS, following the syntax given in the link below. If empty, means the provenance information applies to the whole IDS.
* core_sources/source{:}/global_quantities{:}/electrons/particles: Electron particle source
* core_sources/source{:}/global_quantities{:}/electrons/power: Power coupled to electrons
* core_sources/source{:}/profiles_1d{:}/grid/rho_tor_norm: Normalised toroidal flux coordinate. The normalizing value for rho_tor_norm, is the toroidal flux coordinate at the equilibrium boundary (LCFS or 99.x % of the LCFS in case of a fixed boundary equilibium calculation, see time_slice/boundary/b_flux_pol_norm in the equilibrium IDS) 
* core_sources/source{:}/profiles_1d{:}/grid/rho_tor: Toroidal flux coordinate = sqrt(phi/(pi*b0)), where the toroidal magnetic field, b0, corresponds to that stored in vacuum_toroidal_field/b0 and pi can be found in the IMAS constants
* core_sources/source{:}/profiles_1d{:}/grid/rho_pol_norm: Normalised poloidal flux coordinate = sqrt((psi(rho)-psi(magnetic_axis)) / (psi(LCFS)-psi(magnetic_axis)))
* core_sources/source{:}/profiles_1d{:}/grid/psi: Poloidal magnetic flux
* core_sources/source{:}/profiles_1d{:}/grid/volume: Volume enclosed inside the magnetic surface
* core_sources/source{:}/profiles_1d{:}/grid/area: Cross-sectional area of the flux surface
* core_sources/source{:}/profiles_1d{:}/grid/surface: Surface area of the toroidal flux surface
* core_sources/source{:}/profiles_1d{:}/grid/psi_magnetic_axis: Value of the poloidal magnetic flux at the magnetic axis (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* core_sources/source{:}/profiles_1d{:}/grid/psi_boundary: Value of the poloidal magnetic flux at the plasma boundary (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* core_sources/source{:}/profiles_1d{:}/electrons/particles: Source term for electron density equation
* core_sources/source{:}/profiles_1d{:}/electrons/particles_inside: Electron source inside the flux surface. Cumulative volume integral of the source term for the electron density equation.
* core_sources/source{:}/profiles_1d{:}/electrons/energy: Source term for the electron energy equation
* core_sources/source{:}/profiles_1d{:}/electrons/power_inside: Power coupled to electrons inside the flux surface. Cumulative volume integral of the source term for the electron energy equation
* core_sources/ids_properties/provenance/node{:}/reference{:}/name: Reference name
* core_sources/ids_properties/provenance/node{:}/reference{:}/timestamp: Date and time (UTC) at which the reference was created, expressed in a human readable form (ISO 8601) : the format of the string shall be : YYYY-MM-DDTHH:MM:SSZ. Example : 2020-07-24T14:19:00Z


IDS radiation
--------------

* radiation/time: Generic time
* radiation/ids_properties/comment: Any comment describing the content of this IDS
* radiation/ids_properties/homogeneous_time: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* radiation/ids_properties/source: Source of the data (any comment describing the origin of the data : code, path to diagnostic signals, processing method, ...). Superseeded by the new provenance structure.
* radiation/ids_properties/provider: Name of the person in charge of producing this data
* radiation/ids_properties/creation_date: Date at which this data has been produced
* radiation/vacuum_toroidal_field/r0: Reference major radius where the vacuum toroidal magnetic field is given (usually a fixed position such as the middle of the vessel at the equatorial midplane)
* radiation/vacuum_toroidal_field/b0: Vacuum toroidal field at R0 [T]; Positive sign means anti-clockwise when viewing from above. The product R0B0 must be consistent with the b_tor_vacuum_r field of the tf IDS.
* radiation/code/name: Name of software generating IDS
* radiation/code/commit: Unique commit reference of software
* radiation/code/version: Unique version (tag) of software
* radiation/code/repository: URL of software repository
* radiation/code/parameters: List of the code specific parameters in XML format
* radiation/code/output_flag: Output flag : 0 means the run is successful, other values mean some difficulty has been encountered, the exact meaning is then code specific. Negative values mean the result shall not be used.
* radiation/process{:}/identifier/name: Short string identifier
* radiation/process{:}/identifier/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* radiation/process{:}/identifier/description: Verbose description
* radiation/process{:}/global_quantities{:}/time: Time
* radiation/process{:}/profiles_1d{:}/emissivity_ion_total: Emissivity (summed over ion  species)
* radiation/process{:}/profiles_1d{:}/power_inside_ion_total: Total power from ion species (summed over ion  species) inside the flux surface (volume integral of the emissivity inside the flux surface)
* radiation/process{:}/profiles_1d{:}/emissivity_neutral_total: Emissivity (summed over neutral  species)
* radiation/process{:}/profiles_1d{:}/power_inside_neutral_total: Total power from ion species (summed over neutral species) inside the flux surface (volume integral of the emissivity inside the flux surface)
* radiation/process{:}/profiles_1d{:}/time: Time
* radiation/ids_properties/provenance/node{:}/path: Path of the node within the IDS, following the syntax given in the link below. If empty, means the provenance information applies to the whole IDS.
* radiation/process{:}/global_quantities{:}/inside_lcfs/power: Total power emitted by all species
* radiation/process{:}/global_quantities{:}/inside_lcfs/power_ion_total: Total power emitted by all ion species
* radiation/process{:}/global_quantities{:}/inside_lcfs/power_neutral_total: Total power emitted by all neutral species
* radiation/process{:}/global_quantities{:}/inside_lcfs/power_electrons: Power emitted by electrons
* radiation/process{:}/global_quantities{:}/inside_vessel/power: Total power emitted by all species
* radiation/process{:}/global_quantities{:}/inside_vessel/power_ion_total: Total power emitted by all ion species
* radiation/process{:}/global_quantities{:}/inside_vessel/power_neutral_total: Total power emitted by all neutral species
* radiation/process{:}/global_quantities{:}/inside_vessel/power_electrons: Power emitted by electrons
* radiation/process{:}/profiles_1d{:}/grid/rho_tor_norm: Normalised toroidal flux coordinate. The normalizing value for rho_tor_norm, is the toroidal flux coordinate at the equilibrium boundary (LCFS or 99.x % of the LCFS in case of a fixed boundary equilibium calculation, see time_slice/boundary/b_flux_pol_norm in the equilibrium IDS) 
* radiation/process{:}/profiles_1d{:}/grid/rho_tor: Toroidal flux coordinate = sqrt(phi/(pi*b0)), where the toroidal magnetic field, b0, corresponds to that stored in vacuum_toroidal_field/b0 and pi can be found in the IMAS constants
* radiation/process{:}/profiles_1d{:}/grid/rho_pol_norm: Normalised poloidal flux coordinate = sqrt((psi(rho)-psi(magnetic_axis)) / (psi(LCFS)-psi(magnetic_axis)))
* radiation/process{:}/profiles_1d{:}/grid/psi: Poloidal magnetic flux
* radiation/process{:}/profiles_1d{:}/grid/volume: Volume enclosed inside the magnetic surface
* radiation/process{:}/profiles_1d{:}/grid/area: Cross-sectional area of the flux surface
* radiation/process{:}/profiles_1d{:}/grid/surface: Surface area of the toroidal flux surface
* radiation/process{:}/profiles_1d{:}/grid/psi_magnetic_axis: Value of the poloidal magnetic flux at the magnetic axis (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* radiation/process{:}/profiles_1d{:}/grid/psi_boundary: Value of the poloidal magnetic flux at the plasma boundary (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* radiation/process{:}/profiles_1d{:}/electrons/emissivity: Emissivity from this species
* radiation/process{:}/profiles_1d{:}/electrons/power_inside: Radiated power from inside the flux surface (volume integral of the emissivity inside the flux surface)
* radiation/ids_properties/provenance/node{:}/reference{:}/name: Reference name
* radiation/ids_properties/provenance/node{:}/reference{:}/timestamp: Date and time (UTC) at which the reference was created, expressed in a human readable form (ISO 8601) : the format of the string shall be : YYYY-MM-DDTHH:MM:SSZ. Example : 2020-07-24T14:19:00Z


IDS transport_solver_numerics
------------------------------

* transport_solver_numerics/time: Generic time
* transport_solver_numerics/ids_properties/comment: Any comment describing the content of this IDS
* transport_solver_numerics/ids_properties/homogeneous_time: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* transport_solver_numerics/ids_properties/source: Source of the data (any comment describing the origin of the data : code, path to diagnostic signals, processing method, ...). Superseeded by the new provenance structure.
* transport_solver_numerics/ids_properties/provider: Name of the person in charge of producing this data
* transport_solver_numerics/ids_properties/creation_date: Date at which this data has been produced
* transport_solver_numerics/time_step/data: Data
* transport_solver_numerics/time_step/time: Time
* transport_solver_numerics/solver/name: Short string identifier
* transport_solver_numerics/solver/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* transport_solver_numerics/solver/description: Verbose description
* transport_solver_numerics/primary_coordinate/name: Short string identifier
* transport_solver_numerics/primary_coordinate/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* transport_solver_numerics/primary_coordinate/description: Verbose description
* transport_solver_numerics/solver_1d{:}/drho_tor_dt: Partial derivative of the toroidal flux coordinate profile with respect to time
* transport_solver_numerics/solver_1d{:}/d_dvolume_drho_tor_dt: Partial derivative with respect to time of the derivative of the volume with respect to the toroidal flux coordinate
* transport_solver_numerics/solver_1d{:}/time: Time
* transport_solver_numerics/derivatives_1d{:}/dpsi_dt: Derivative of the poloidal flux profile with respect to time
* transport_solver_numerics/derivatives_1d{:}/dpsi_dt_cphi: Derivative of the poloidal flux profile with respect to time, at constant toroidal flux
* transport_solver_numerics/derivatives_1d{:}/dpsi_dt_crho_tor_norm: Derivative of the poloidal flux profile with respect to time, at constant normalised toroidal flux coordinate
* transport_solver_numerics/derivatives_1d{:}/drho_tor_dt: Partial derivative of the toroidal flux coordinate profile with respect to time
* transport_solver_numerics/derivatives_1d{:}/d_dvolume_drho_tor_dt: Partial derivative with respect to time of the derivative of the volume with respect to the toroidal flux coordinate
* transport_solver_numerics/derivatives_1d{:}/dpsi_drho_tor: Derivative of the poloidal flux profile with respect to the toroidal flux coordinate
* transport_solver_numerics/derivatives_1d{:}/d2psi_drho_tor2: Second derivative of the poloidal flux profile with respect to the toroidal flux coordinate
* transport_solver_numerics/derivatives_1d{:}/time: Time
* transport_solver_numerics/boundary_conditions_1d{:}/time: Time
* transport_solver_numerics/vacuum_toroidal_field/r0: Reference major radius where the vacuum toroidal magnetic field is given (usually a fixed position such as the middle of the vessel at the equatorial midplane)
* transport_solver_numerics/vacuum_toroidal_field/b0: Vacuum toroidal field at R0 [T]; Positive sign means anti-clockwise when viewing from above. The product R0B0 must be consistent with the b_tor_vacuum_r field of the tf IDS.
* transport_solver_numerics/restart_files{:}/time: Time
* transport_solver_numerics/code/name: Name of software generating IDS
* transport_solver_numerics/code/commit: Unique commit reference of software
* transport_solver_numerics/code/version: Unique version (tag) of software
* transport_solver_numerics/code/repository: URL of software repository
* transport_solver_numerics/code/parameters: List of the code specific parameters in XML format
* transport_solver_numerics/code/output_flag: Output flag : 0 means the run is successful, other values mean some difficulty has been encountered, the exact meaning is then code specific. Negative values mean the result shall not be used.
* transport_solver_numerics/solver_1d{:}/grid/rho_tor_norm: Normalised toroidal flux coordinate. The normalizing value for rho_tor_norm, is the toroidal flux coordinate at the equilibrium boundary (LCFS or 99.x % of the LCFS in case of a fixed boundary equilibium calculation, see time_slice/boundary/b_flux_pol_norm in the equilibrium IDS) 
* transport_solver_numerics/solver_1d{:}/grid/rho_tor: Toroidal flux coordinate = sqrt(phi/(pi*b0)), where the toroidal magnetic field, b0, corresponds to that stored in vacuum_toroidal_field/b0 and pi can be found in the IMAS constants
* transport_solver_numerics/solver_1d{:}/grid/rho_pol_norm: Normalised poloidal flux coordinate = sqrt((psi(rho)-psi(magnetic_axis)) / (psi(LCFS)-psi(magnetic_axis)))
* transport_solver_numerics/solver_1d{:}/grid/psi: Poloidal magnetic flux
* transport_solver_numerics/solver_1d{:}/grid/volume: Volume enclosed inside the magnetic surface
* transport_solver_numerics/solver_1d{:}/grid/area: Cross-sectional area of the flux surface
* transport_solver_numerics/solver_1d{:}/grid/surface: Surface area of the toroidal flux surface
* transport_solver_numerics/solver_1d{:}/grid/psi_magnetic_axis: Value of the poloidal magnetic flux at the magnetic axis (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* transport_solver_numerics/solver_1d{:}/grid/psi_boundary: Value of the poloidal magnetic flux at the plasma boundary (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* transport_solver_numerics/derivatives_1d{:}/grid/rho_tor_norm: Normalised toroidal flux coordinate. The normalizing value for rho_tor_norm, is the toroidal flux coordinate at the equilibrium boundary (LCFS or 99.x % of the LCFS in case of a fixed boundary equilibium calculation, see time_slice/boundary/b_flux_pol_norm in the equilibrium IDS) 
* transport_solver_numerics/derivatives_1d{:}/grid/rho_tor: Toroidal flux coordinate = sqrt(phi/(pi*b0)), where the toroidal magnetic field, b0, corresponds to that stored in vacuum_toroidal_field/b0 and pi can be found in the IMAS constants
* transport_solver_numerics/derivatives_1d{:}/grid/rho_pol_norm: Normalised poloidal flux coordinate = sqrt((psi(rho)-psi(magnetic_axis)) / (psi(LCFS)-psi(magnetic_axis)))
* transport_solver_numerics/derivatives_1d{:}/grid/psi: Poloidal magnetic flux
* transport_solver_numerics/derivatives_1d{:}/grid/volume: Volume enclosed inside the magnetic surface
* transport_solver_numerics/derivatives_1d{:}/grid/area: Cross-sectional area of the flux surface
* transport_solver_numerics/derivatives_1d{:}/grid/surface: Surface area of the toroidal flux surface
* transport_solver_numerics/derivatives_1d{:}/grid/psi_magnetic_axis: Value of the poloidal magnetic flux at the magnetic axis (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* transport_solver_numerics/derivatives_1d{:}/grid/psi_boundary: Value of the poloidal magnetic flux at the plasma boundary (useful to normalize the psi array values when the radial grid doesn't go from the magnetic axis to the plasma boundary)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/a: Mass of atom
* transport_solver_numerics/derivatives_1d{:}/ion{:}/z_ion: Ion charge (of the dominant ionisation state; lumped ions are allowed)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/z_n: Nuclear charge
* transport_solver_numerics/derivatives_1d{:}/ion{:}/label: String identifying ion (e.g. H+, D+, T+, He+2, C+, ...)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/name: String identifying ion (e.g. H+, D+, T+, He+2, C+, ...)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/multiple_states_flag: Multiple state calculation flag : 0-Only one state is considered; 1-Multiple states are considered and are described in the state structure
* transport_solver_numerics/derivatives_1d{:}/d_drho_tor_norm/n_i_total_over_n_e: Ratio of total ion density (sum over species and charge states) over electron density. (thermal+non-thermal)
* transport_solver_numerics/derivatives_1d{:}/d_drho_tor_norm/pressure_ion_total: Total thermal ion pressure
* transport_solver_numerics/derivatives_1d{:}/d2_drho_tor_norm2/n_i_total_over_n_e: Ratio of total ion density (sum over species and charge states) over electron density. (thermal+non-thermal)
* transport_solver_numerics/derivatives_1d{:}/d2_drho_tor_norm2/pressure_ion_total: Total thermal ion pressure
* transport_solver_numerics/derivatives_1d{:}/d_dt/n_i_total_over_n_e: Ratio of total ion density (sum over species and charge states) over electron density. (thermal+non-thermal)
* transport_solver_numerics/derivatives_1d{:}/d_dt/pressure_ion_total: Total thermal ion pressure
* transport_solver_numerics/boundary_conditions_1d{:}/current/value: Value of the boundary condition. For ID = 1 to 3, only the first position in the vector is used. For ID = 5, all three positions are used, meaning respectively a1, a2, a3. 
* transport_solver_numerics/boundary_conditions_1d{:}/current/rho_tor_norm: Position, in normalised toroidal flux, at which the boundary condition is imposed. Outside this position, the value of the data are considered to be prescribed.
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/a: Mass of atom
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/z_ion: Ion charge (of the dominant ionisation state; lumped ions are allowed)
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/z_n: Nuclear charge
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/label: String identifying ion (e.g. H+, D+, T+, He+2, C+, ...)
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/name: String identifying ion (e.g. H+, D+, T+, He+2, C+, ...)
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/multiple_states_flag: Multiple states calculation flag : 0-Only one state is considered; 1-Multiple states are considered and are described in the state structure
* transport_solver_numerics/boundary_conditions_1d{:}/energy_ion_total/value: Value of the boundary condition. For ID = 1 to 4, only the first position in the vector is used. For ID = 5, all three positions are used, meaning respectively a1, a2, a3. 
* transport_solver_numerics/boundary_conditions_1d{:}/energy_ion_total/rho_tor_norm: Position, in normalised toroidal flux, at which the boundary condition is imposed. Outside this position, the value of the data are considered to be prescribed.
* transport_solver_numerics/boundary_conditions_1d{:}/momentum_tor/value: Value of the boundary condition. For ID = 1 to 4, only the first position in the vector is used. For ID = 5, all three positions are used, meaning respectively a1, a2, a3. 
* transport_solver_numerics/boundary_conditions_1d{:}/momentum_tor/rho_tor_norm: Position, in normalised toroidal flux, at which the boundary condition is imposed. Outside this position, the value of the data are considered to be prescribed.
* transport_solver_numerics/convergence/time_step/data: Data
* transport_solver_numerics/convergence/time_step/time: Time
* transport_solver_numerics/convergence/equations{:}/time: Time
* transport_solver_numerics/ids_properties/provenance/node{:}/path: Path of the node within the IDS, following the syntax given in the link below. If empty, means the provenance information applies to the whole IDS.
* transport_solver_numerics/solver_1d{:}/equation{:}/primary_quantity/profile: Profile of the primary quantity
* transport_solver_numerics/solver_1d{:}/equation{:}/primary_quantity/d_dr: Radial derivative with respect to the primary coordinate
* transport_solver_numerics/solver_1d{:}/equation{:}/primary_quantity/d2_dr2: Second order radial derivative with respect to the primary coordinate
* transport_solver_numerics/solver_1d{:}/equation{:}/primary_quantity/d_dt: Time derivative
* transport_solver_numerics/solver_1d{:}/equation{:}/primary_quantity/d_dt_cphi: Derivative with respect to time, at constant toroidal flux (for current diffusion equation)
* transport_solver_numerics/solver_1d{:}/equation{:}/primary_quantity/d_dt_cr: Derivative with respect to time, at constant primary coordinate coordinate (for current diffusion equation)
* transport_solver_numerics/solver_1d{:}/equation{:}/boundary_condition{:}/value: Value of the boundary condition. For type/index = 1 to 3, only the first position in the vector is used. For type/index = 5, all three positions are used, meaning respectively a1, a2, a3. 
* transport_solver_numerics/solver_1d{:}/equation{:}/boundary_condition{:}/position: Position, in terms of the primary coordinate, at which the boundary condition is imposed. Outside this position, the value of the data are considered to be prescribed (in case of a single boundary condition).
* transport_solver_numerics/solver_1d{:}/equation{:}/coefficient{:}/profile: Radial profile of the numerical coefficient
* transport_solver_numerics/solver_1d{:}/equation{:}/convergence/iterations_n: Number of iterations carried out in the convergence loop
* transport_solver_numerics/derivatives_1d{:}/electrons/d_drho_tor_norm/temperature: Temperature
* transport_solver_numerics/derivatives_1d{:}/electrons/d_drho_tor_norm/density: Density (thermal+non-thermal)
* transport_solver_numerics/derivatives_1d{:}/electrons/d_drho_tor_norm/density_fast: Density of fast (non-thermal) particles
* transport_solver_numerics/derivatives_1d{:}/electrons/d_drho_tor_norm/pressure: Pressure
* transport_solver_numerics/derivatives_1d{:}/electrons/d_drho_tor_norm/pressure_fast_perpendicular: Fast (non-thermal) perpendicular pressure
* transport_solver_numerics/derivatives_1d{:}/electrons/d_drho_tor_norm/pressure_fast_parallel: Fast (non-thermal) parallel pressure
* transport_solver_numerics/derivatives_1d{:}/electrons/d2_drho_tor_norm2/temperature: Temperature
* transport_solver_numerics/derivatives_1d{:}/electrons/d2_drho_tor_norm2/density: Density (thermal+non-thermal)
* transport_solver_numerics/derivatives_1d{:}/electrons/d2_drho_tor_norm2/density_fast: Density of fast (non-thermal) particles
* transport_solver_numerics/derivatives_1d{:}/electrons/d2_drho_tor_norm2/pressure: Pressure
* transport_solver_numerics/derivatives_1d{:}/electrons/d2_drho_tor_norm2/pressure_fast_perpendicular: Fast (non-thermal) perpendicular pressure
* transport_solver_numerics/derivatives_1d{:}/electrons/d2_drho_tor_norm2/pressure_fast_parallel: Fast (non-thermal) parallel pressure
* transport_solver_numerics/derivatives_1d{:}/electrons/d_dt/temperature: Temperature
* transport_solver_numerics/derivatives_1d{:}/electrons/d_dt/density: Density (thermal+non-thermal)
* transport_solver_numerics/derivatives_1d{:}/electrons/d_dt/density_fast: Density of fast (non-thermal) particles
* transport_solver_numerics/derivatives_1d{:}/electrons/d_dt/pressure: Pressure
* transport_solver_numerics/derivatives_1d{:}/electrons/d_dt/pressure_fast_perpendicular: Fast (non-thermal) perpendicular pressure
* transport_solver_numerics/derivatives_1d{:}/electrons/d_dt/pressure_fast_parallel: Fast (non-thermal) parallel pressure
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_drho_tor_norm/temperature: Temperature (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_drho_tor_norm/density: Density (thermal+non-thermal) (sum over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_drho_tor_norm/density_fast: Density of fast (non-thermal) particles (sum over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_drho_tor_norm/pressure: Pressure (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_drho_tor_norm/pressure_fast_perpendicular: Fast (non-thermal) perpendicular pressure  (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_drho_tor_norm/pressure_fast_parallel: Fast (non-thermal) parallel pressure  (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_drho_tor_norm/velocity_tor: Toroidal velocity (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_drho_tor_norm/velocity_pol: Poloidal velocity (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d2_drho_tor_norm2/temperature: Temperature (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d2_drho_tor_norm2/density: Density (thermal+non-thermal) (sum over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d2_drho_tor_norm2/density_fast: Density of fast (non-thermal) particles (sum over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d2_drho_tor_norm2/pressure: Pressure (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d2_drho_tor_norm2/pressure_fast_perpendicular: Fast (non-thermal) perpendicular pressure  (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d2_drho_tor_norm2/pressure_fast_parallel: Fast (non-thermal) parallel pressure  (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d2_drho_tor_norm2/velocity_tor: Toroidal velocity (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d2_drho_tor_norm2/velocity_pol: Poloidal velocity (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_dt/temperature: Temperature (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_dt/density: Density (thermal+non-thermal) (sum over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_dt/density_fast: Density of fast (non-thermal) particles (sum over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_dt/pressure: Pressure (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_dt/pressure_fast_perpendicular: Fast (non-thermal) perpendicular pressure  (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_dt/pressure_fast_parallel: Fast (non-thermal) parallel pressure  (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_dt/velocity_tor: Toroidal velocity (average over charge states when multiple charge states are considered)
* transport_solver_numerics/derivatives_1d{:}/ion{:}/d_dt/velocity_pol: Poloidal velocity (average over charge states when multiple charge states are considered)
* transport_solver_numerics/boundary_conditions_1d{:}/current/identifier/name: Short string identifier
* transport_solver_numerics/boundary_conditions_1d{:}/current/identifier/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* transport_solver_numerics/boundary_conditions_1d{:}/current/identifier/description: Verbose description
* transport_solver_numerics/boundary_conditions_1d{:}/electrons/particles/value: Value of the boundary condition. For ID = 1 to 4, only the first position in the vector is used. For ID = 5, all three positions are used, meaning respectively a1, a2, a3. 
* transport_solver_numerics/boundary_conditions_1d{:}/electrons/particles/rho_tor_norm: Position, in normalised toroidal flux, at which the boundary condition is imposed. Outside this position, the value of the data are considered to be prescribed.
* transport_solver_numerics/boundary_conditions_1d{:}/electrons/energy/value: Value of the boundary condition. For ID = 1 to 4, only the first position in the vector is used. For ID = 5, all three positions are used, meaning respectively a1, a2, a3. 
* transport_solver_numerics/boundary_conditions_1d{:}/electrons/energy/rho_tor_norm: Position, in normalised toroidal flux, at which the boundary condition is imposed. Outside this position, the value of the data are considered to be prescribed.
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/particles/value: Value of the boundary condition. For ID = 1 to 4, only the first position in the vector is used. For ID = 5, all three positions are used, meaning respectively a1, a2, a3. 
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/particles/rho_tor_norm: Position, in normalised toroidal flux, at which the boundary condition is imposed. Outside this position, the value of the data are considered to be prescribed.
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/energy/value: Value of the boundary condition. For ID = 1 to 4, only the first position in the vector is used. For ID = 5, all three positions are used, meaning respectively a1, a2, a3. 
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/energy/rho_tor_norm: Position, in normalised toroidal flux, at which the boundary condition is imposed. Outside this position, the value of the data are considered to be prescribed.
* transport_solver_numerics/boundary_conditions_1d{:}/energy_ion_total/identifier/name: Short string identifier
* transport_solver_numerics/boundary_conditions_1d{:}/energy_ion_total/identifier/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* transport_solver_numerics/boundary_conditions_1d{:}/energy_ion_total/identifier/description: Verbose description
* transport_solver_numerics/boundary_conditions_1d{:}/momentum_tor/identifier/name: Short string identifier
* transport_solver_numerics/boundary_conditions_1d{:}/momentum_tor/identifier/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* transport_solver_numerics/boundary_conditions_1d{:}/momentum_tor/identifier/description: Verbose description
* transport_solver_numerics/convergence/equations{:}/current/iterations_n: Number of iterations carried out in the convergence loop
* transport_solver_numerics/convergence/equations{:}/ion{:}/a: Mass of atom
* transport_solver_numerics/convergence/equations{:}/ion{:}/z_ion: Ion charge (of the dominant ionisation state; lumped ions are allowed)
* transport_solver_numerics/convergence/equations{:}/ion{:}/z_n: Nuclear charge
* transport_solver_numerics/convergence/equations{:}/ion{:}/label: String identifying ion (e.g. H+, D+, T+, He+2, C+, ...)
* transport_solver_numerics/convergence/equations{:}/ion{:}/name: String identifying ion (e.g. H+, D+, T+, He+2, C+, ...)
* transport_solver_numerics/convergence/equations{:}/ion{:}/multiple_states_flag: Multiple state calculation flag : 0-Only one state is considered; 1-Multiple states are considered and are described in the state structure
* transport_solver_numerics/convergence/equations{:}/energy_ion_total/iterations_n: Number of iterations carried out in the convergence loop
* transport_solver_numerics/ids_properties/provenance/node{:}/reference{:}/name: Reference name
* transport_solver_numerics/ids_properties/provenance/node{:}/reference{:}/timestamp: Date and time (UTC) at which the reference was created, expressed in a human readable form (ISO 8601) : the format of the string shall be : YYYY-MM-DDTHH:MM:SSZ. Example : 2020-07-24T14:19:00Z
* transport_solver_numerics/solver_1d{:}/equation{:}/primary_quantity/identifier/name: Short string identifier
* transport_solver_numerics/solver_1d{:}/equation{:}/primary_quantity/identifier/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* transport_solver_numerics/solver_1d{:}/equation{:}/primary_quantity/identifier/description: Verbose description
* transport_solver_numerics/solver_1d{:}/equation{:}/boundary_condition{:}/type/name: Short string identifier
* transport_solver_numerics/solver_1d{:}/equation{:}/boundary_condition{:}/type/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* transport_solver_numerics/solver_1d{:}/equation{:}/boundary_condition{:}/type/description: Verbose description
* transport_solver_numerics/boundary_conditions_1d{:}/electrons/particles/identifier/name: Short string identifier
* transport_solver_numerics/boundary_conditions_1d{:}/electrons/particles/identifier/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* transport_solver_numerics/boundary_conditions_1d{:}/electrons/particles/identifier/description: Verbose description
* transport_solver_numerics/boundary_conditions_1d{:}/electrons/energy/identifier/name: Short string identifier
* transport_solver_numerics/boundary_conditions_1d{:}/electrons/energy/identifier/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* transport_solver_numerics/boundary_conditions_1d{:}/electrons/energy/identifier/description: Verbose description
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/particles/identifier/name: Short string identifier
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/particles/identifier/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/particles/identifier/description: Verbose description
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/energy/identifier/name: Short string identifier
* transport_solver_numerics/boundary_conditions_1d{:}/ion{:}/energy/identifier/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* transport_solver_numerics/convergence/equations{:}/electrons/particles/iterations_n: Number of iterations carried out in the convergence loop
* transport_solver_numerics/convergence/equations{:}/electrons/energy/iterations_n: Number of iterations carried out in the convergence loop
* transport_solver_numerics/convergence/equations{:}/ion{:}/particles/iterations_n: Number of iterations carried out in the convergence loop
* transport_solver_numerics/convergence/equations{:}/ion{:}/energy/iterations_n: Number of iterations carried out in the convergence loop


IDS equilibrium
----------------

* equilibrium/time: Generic time
* equilibrium/ids_properties/comment: Any comment describing the content of this IDS
* equilibrium/ids_properties/homogeneous_time: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* equilibrium/ids_properties/source: Source of the data (any comment describing the origin of the data : code, path to diagnostic signals, processing method, ...). Superseeded by the new provenance structure.
* equilibrium/ids_properties/provider: Name of the person in charge of producing this data
* equilibrium/ids_properties/creation_date: Date at which this data has been produced
* equilibrium/vacuum_toroidal_field/r0: Reference major radius where the vacuum toroidal magnetic field is given (usually a fixed position such as the middle of the vessel at the equatorial midplane)
* equilibrium/vacuum_toroidal_field/b0: Vacuum toroidal field at R0 [T]; Positive sign means anti-clockwise when viewing from above. The product R0B0 must be consistent with the b_tor_vacuum_r field of the tf IDS.
* equilibrium/grids_ggd{:}/time: Time
* equilibrium/time_slice{:}/time: Time
* equilibrium/code/name: Name of software generating IDS
* equilibrium/code/commit: Unique commit reference of software
* equilibrium/code/version: Unique version (tag) of software
* equilibrium/code/repository: URL of software repository
* equilibrium/code/parameters: List of the code specific parameters in XML format
* equilibrium/code/output_flag: Output flag : 0 means the run is successful, other values mean some difficulty has been encountered, the exact meaning is then code specific. Negative values mean the result shall not be used.
* equilibrium/time_slice{:}/boundary/type: 0 (limiter) or 1 (diverted)
* equilibrium/time_slice{:}/boundary/psi: Value of the poloidal flux at which the boundary is taken
* equilibrium/time_slice{:}/boundary/minor_radius: Minor radius of the plasma boundary (defined as (Rmax-Rmin) / 2 of the boundary)
* equilibrium/time_slice{:}/boundary/elongation: Elongation of the plasma boundary
* equilibrium/time_slice{:}/boundary/elongation_upper: Elongation (upper half w.r.t. geometric axis) of the plasma boundary
* equilibrium/time_slice{:}/boundary/elongation_lower: Elongation (lower half w.r.t. geometric axis) of the plasma boundary
* equilibrium/time_slice{:}/boundary/triangularity_upper: Upper triangularity of the plasma boundary
* equilibrium/time_slice{:}/boundary/triangularity_lower: Lower triangularity of the plasma boundary
* equilibrium/time_slice{:}/boundary/squareness_upper_inner: Upper inner squareness of the plasma boundary (definition from T. Luce, Plasma Phys. Control. Fusion 55 (2013) 095009)
* equilibrium/time_slice{:}/boundary/squareness_upper_outer: Upper outer squareness of the plasma boundary (definition from T. Luce, Plasma Phys. Control. Fusion 55 (2013) 095009)
* equilibrium/time_slice{:}/boundary/squareness_lower_inner: Lower inner squareness of the plasma boundary (definition from T. Luce, Plasma Phys. Control. Fusion 55 (2013) 095009)
* equilibrium/time_slice{:}/boundary/squareness_lower_outer: Lower outer squareness of the plasma boundary (definition from T. Luce, Plasma Phys. Control. Fusion 55 (2013) 095009)
* equilibrium/time_slice{:}/boundary_separatrix/type: 0 (limiter) or 1 (diverted)
* equilibrium/time_slice{:}/boundary_separatrix/psi: Value of the poloidal flux at the separatrix
* equilibrium/time_slice{:}/boundary_separatrix/minor_radius: Minor radius of the plasma boundary (defined as (Rmax-Rmin) / 2 of the boundary)
* equilibrium/time_slice{:}/boundary_separatrix/elongation: Elongation of the plasma boundary
* equilibrium/time_slice{:}/boundary_separatrix/elongation_upper: Elongation (upper half w.r.t. geometric axis) of the plasma boundary
* equilibrium/time_slice{:}/boundary_separatrix/elongation_lower: Elongation (lower half w.r.t. geometric axis) of the plasma boundary
* equilibrium/time_slice{:}/boundary_separatrix/triangularity_upper: Upper triangularity of the plasma boundary
* equilibrium/time_slice{:}/boundary_separatrix/triangularity_lower: Lower triangularity of the plasma boundary
* equilibrium/time_slice{:}/boundary_separatrix/squareness_upper_inner: Upper inner squareness of the plasma boundary (definition from T. Luce, Plasma Phys. Control. Fusion 55 (2013) 095009)
* equilibrium/time_slice{:}/boundary_separatrix/squareness_upper_outer: Upper outer squareness of the plasma boundary (definition from T. Luce, Plasma Phys. Control. Fusion 55 (2013) 095009)
* equilibrium/time_slice{:}/boundary_separatrix/squareness_lower_inner: Lower inner squareness of the plasma boundary (definition from T. Luce, Plasma Phys. Control. Fusion 55 (2013) 095009)
* equilibrium/time_slice{:}/boundary_separatrix/squareness_lower_outer: Lower outer squareness of the plasma boundary (definition from T. Luce, Plasma Phys. Control. Fusion 55 (2013) 095009)
* equilibrium/time_slice{:}/global_quantities/beta_pol: Poloidal beta. Defined as betap = 4 int(p dV) / [R_0 * mu_0 * Ip^2]
* equilibrium/time_slice{:}/global_quantities/beta_tor: Toroidal beta, defined as the volume-averaged total perpendicular pressure divided by (B0^2/(2*mu0)), i.e. beta_toroidal = 2 mu0 int(p dV) / V / B0^2
* equilibrium/time_slice{:}/global_quantities/beta_normal: Normalised toroidal beta, defined as 100 * beta_tor * a[m] * B0 [T] / ip [MA] 
* equilibrium/time_slice{:}/global_quantities/ip: Plasma current (toroidal component). Positive sign means anti-clockwise when viewed from above.
* equilibrium/time_slice{:}/global_quantities/ip_error_upper: Upper error for "ip"
* equilibrium/time_slice{:}/global_quantities/ip_error_lower: Lower error for "ip"
* equilibrium/time_slice{:}/global_quantities/li_3: Internal inductance
* equilibrium/time_slice{:}/global_quantities/volume: Total plasma volume
* equilibrium/time_slice{:}/global_quantities/area: Area of the LCFS poloidal cross section
* equilibrium/time_slice{:}/global_quantities/surface: Surface area of the toroidal flux surface
* equilibrium/time_slice{:}/global_quantities/length_pol: Poloidal length of the magnetic surface
* equilibrium/time_slice{:}/global_quantities/psi_axis: Poloidal flux at the magnetic axis
* equilibrium/time_slice{:}/global_quantities/psi_boundary: Poloidal flux at the selected plasma boundary 
* equilibrium/time_slice{:}/global_quantities/psi_boundary_error_upper: Upper error for "psi_boundary"
* equilibrium/time_slice{:}/global_quantities/psi_boundary_error_lower: Lower error for "psi_boundary"
* equilibrium/time_slice{:}/global_quantities/q_axis: q at the magnetic axis
* equilibrium/time_slice{:}/global_quantities/q_95: q at the 95% poloidal flux surface (IMAS uses COCOS=11: only positive when toroidal current and magnetic field are in same direction)
* equilibrium/time_slice{:}/global_quantities/energy_mhd: Plasma energy content = 3/2 * int(p,dV) with p being the total pressure (thermal + fast particles) [J]. Time-dependent; Scalar
* equilibrium/time_slice{:}/global_quantities/w_mhd: Plasma energy content = 3/2 * int(p,dV) with p being the total pressure (thermal + fast particles) [J]. Time-dependent; Scalar
* equilibrium/time_slice{:}/global_quantities/plasma_resistance: Plasma resistance = int(e_field.j.dV) / Ip^2
* equilibrium/time_slice{:}/profiles_1d/psi: Poloidal flux
* equilibrium/time_slice{:}/profiles_1d/phi: Toroidal flux
* equilibrium/time_slice{:}/profiles_1d/pressure: Pressure
* equilibrium/time_slice{:}/profiles_1d/f: Diamagnetic function (F=R B_Phi)
* equilibrium/time_slice{:}/profiles_1d/dpressure_dpsi: Derivative of pressure w.r.t. psi
* equilibrium/time_slice{:}/profiles_1d/f_df_dpsi: Derivative of F w.r.t. Psi, multiplied with F
* equilibrium/time_slice{:}/profiles_1d/j_tor: Flux surface averaged toroidal current density = average(j_tor/R) / average(1/R)
* equilibrium/time_slice{:}/profiles_1d/j_phi: Flux surface averaged toroidal current density = average(j_tor/R) / average(1/R)
* equilibrium/time_slice{:}/profiles_1d/j_parallel: Flux surface averaged approximation to parallel current density = average(j.B) / B0, where B0 = /vacuum_toroidal_field/b0
* equilibrium/time_slice{:}/profiles_1d/q: Safety factor (IMAS uses COCOS=11: only positive when toroidal current and magnetic field are in same direction)
* equilibrium/time_slice{:}/profiles_1d/magnetic_shear: Magnetic shear, defined as rho_tor/q . dq/drho_tor
* equilibrium/time_slice{:}/profiles_1d/r_inboard: Radial coordinate (major radius) on the inboard side of the magnetic axis
* equilibrium/time_slice{:}/profiles_1d/r_outboard: Radial coordinate (major radius) on the outboard side of the magnetic axis
* equilibrium/time_slice{:}/profiles_1d/rho_tor: Toroidal flux coordinate = sqrt(phi/(pi*b0)), where the toroidal flux, phi, corresponds to time_slice/profiles_1d/phi, the toroidal magnetic field, b0, corresponds to that stored in vacuum_toroidal_field/b0 and pi can be found in the IMAS constants
* equilibrium/time_slice{:}/profiles_1d/rho_tor_norm: Normalised toroidal flux coordinate. The normalizing value for rho_tor_norm, is the toroidal flux coordinate at the equilibrium boundary (LCFS or 99.x % of the LCFS in case of a fixed boundary equilibium calculation) 
* equilibrium/time_slice{:}/profiles_1d/dpsi_drho_tor: Derivative of Psi with respect to Rho_Tor
* equilibrium/time_slice{:}/profiles_1d/elongation: Elongation
* equilibrium/time_slice{:}/profiles_1d/triangularity_upper: Upper triangularity w.r.t. magnetic axis
* equilibrium/time_slice{:}/profiles_1d/triangularity_lower: Lower triangularity w.r.t. magnetic axis
* equilibrium/time_slice{:}/profiles_1d/volume: Volume enclosed in the flux surface
* equilibrium/time_slice{:}/profiles_1d/dvolume_dpsi: Radial derivative of the volume enclosed in the flux surface with respect to Psi
* equilibrium/time_slice{:}/profiles_1d/dvolume_drho_tor: Radial derivative of the volume enclosed in the flux surface with respect to Rho_Tor
* equilibrium/time_slice{:}/profiles_1d/area: Cross-sectional area of the flux surface
* equilibrium/time_slice{:}/profiles_1d/darea_dpsi: Radial derivative of the cross-sectional area of the flux surface with respect to psi
* equilibrium/time_slice{:}/profiles_1d/darea_drho_tor: Radial derivative of the cross-sectional area of the flux surface with respect to rho_tor
* equilibrium/time_slice{:}/profiles_1d/surface: Surface area of the toroidal flux surface
* equilibrium/time_slice{:}/profiles_1d/trapped_fraction: Trapped particle fraction
* equilibrium/time_slice{:}/profiles_1d/gm1: Flux surface averaged 1/R^2
* equilibrium/time_slice{:}/profiles_1d/gm2: Flux surface averaged |grad_rho_tor|^2/R^2
* equilibrium/time_slice{:}/profiles_1d/gm3: Flux surface averaged |grad_rho_tor|^2
* equilibrium/time_slice{:}/profiles_1d/gm4: Flux surface averaged 1/B^2
* equilibrium/time_slice{:}/profiles_1d/gm5: Flux surface averaged B^2
* equilibrium/time_slice{:}/profiles_1d/gm6: Flux surface averaged |grad_rho_tor|^2/B^2
* equilibrium/time_slice{:}/profiles_1d/gm7: Flux surface averaged |grad_rho_tor|
* equilibrium/time_slice{:}/profiles_1d/gm8: Flux surface averaged R
* equilibrium/time_slice{:}/profiles_1d/gm9: Flux surface averaged 1/R
* equilibrium/time_slice{:}/profiles_1d/b_average: Flux surface averaged B
* equilibrium/time_slice{:}/profiles_1d/b_field_average: Flux surface averaged modulus of B (always positive, irrespective of the sign convention for the B-field direction).
* equilibrium/time_slice{:}/profiles_1d/b_min: Minimum(B) on the flux surface
* equilibrium/time_slice{:}/profiles_1d/b_field_min: Minimum(modulus(B)) on the flux surface (always positive, irrespective of the sign convention for the B-field direction)
* equilibrium/time_slice{:}/profiles_1d/b_max: Maximum(B) on the flux surface
* equilibrium/time_slice{:}/profiles_1d/b_field_max: Maximum(modulus(B)) on the flux surface (always positive, irrespective of the sign convention for the B-field direction)
* equilibrium/time_slice{:}/coordinate_system/r: Values of the major radius on the grid
* equilibrium/time_slice{:}/coordinate_system/z: Values of the Height on the grid
* equilibrium/time_slice{:}/coordinate_system/jacobian: Absolute value of the jacobian of the coordinate system
* equilibrium/time_slice{:}/coordinate_system/tensor_covariant: Covariant metric tensor on every point of the grid described by grid_type
* equilibrium/time_slice{:}/coordinate_system/tensor_contravariant: Contravariant metric tensor on every point of the grid described by grid_type
* equilibrium/time_slice{:}/coordinate_system/g11_covariant: metric coefficients g11,  covariant metric tensor for the grid described by grid_type
* equilibrium/time_slice{:}/coordinate_system/g12_covariant: metric coefficients g12, covariant metric tensor for the grid described by grid_type
* equilibrium/time_slice{:}/coordinate_system/g13_covariant: metric coefficients g13, covariant metric tensor for the grid described by grid_type
* equilibrium/time_slice{:}/coordinate_system/g22_covariant: metric coefficients g22, covariant metric tensor for the grid described by grid_type
* equilibrium/time_slice{:}/coordinate_system/g23_covariant: metric coefficients g23,  covariant metric tensor for the grid described by grid_type
* equilibrium/time_slice{:}/coordinate_system/g33_covariant: metric coefficients g33, covariant metric tensor for the grid described by grid_type
* equilibrium/time_slice{:}/coordinate_system/g11_contravariant: metric coefficients g11, contravariant metric tensor for the grid described by grid_type
* equilibrium/time_slice{:}/coordinate_system/g12_contravariant: metric coefficients g12, contravariant metric tensor for the grid described by grid_type
* equilibrium/time_slice{:}/coordinate_system/g13_contravariant: metric coefficients g13,  contravariant metric tensor for the grid described by grid_type
* equilibrium/time_slice{:}/coordinate_system/g22_contravariant: metric coefficients g22,  contravariant metric tensor for the grid described by grid_type
* equilibrium/time_slice{:}/coordinate_system/g23_contravariant: metric coefficients g23, contravariant metric tensor for the grid described by grid_type
* equilibrium/time_slice{:}/coordinate_system/g33_contravariant: metric coefficients g33,  contravariant metric tensor for the grid described by grid_type
* equilibrium/time_slice{:}/profiles_2d{:}/r: Values of the major radius on the grid
* equilibrium/time_slice{:}/profiles_2d{:}/z: Values of the Height on the grid
* equilibrium/time_slice{:}/profiles_2d{:}/psi: Values of the poloidal flux at the grid in the poloidal plane
* equilibrium/time_slice{:}/profiles_2d{:}/theta: Values of the poloidal angle on the grid
* equilibrium/time_slice{:}/profiles_2d{:}/phi: Toroidal flux
* equilibrium/time_slice{:}/profiles_2d{:}/j_tor: Toroidal plasma current density
* equilibrium/time_slice{:}/profiles_2d{:}/j_phi: Toroidal plasma current density
* equilibrium/time_slice{:}/profiles_2d{:}/j_parallel: Defined as (j.B)/B0 where j and B are the current density and magnetic field vectors and B0 is the (signed) vacuum toroidal magnetic field strength at the geometric reference point (R0,Z0). It is formally not the component of the plasma current density parallel to the magnetic field
* equilibrium/time_slice{:}/profiles_2d{:}/b_r: R component of the poloidal magnetic field
* equilibrium/time_slice{:}/profiles_2d{:}/b_field_r: R component of the poloidal magnetic field
* equilibrium/time_slice{:}/profiles_2d{:}/b_field_phi: Toroidal component of the magnetic field
* equilibrium/time_slice{:}/profiles_2d{:}/b_z: Z component of the poloidal magnetic field
* equilibrium/time_slice{:}/profiles_2d{:}/b_field_z: Z component of the poloidal magnetic field
* equilibrium/time_slice{:}/profiles_2d{:}/b_tor: Toroidal component of the magnetic field
* equilibrium/time_slice{:}/profiles_2d{:}/b_field_tor: Toroidal component of the magnetic field
* equilibrium/ids_properties/provenance/node{:}/path: Path of the node within the IDS, following the syntax given in the link below. If empty, means the provenance information applies to the whole IDS.
* equilibrium/grids_ggd{:}/grid{:}/identifier/name: Short string identifier
* equilibrium/grids_ggd{:}/grid{:}/identifier/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* equilibrium/grids_ggd{:}/grid{:}/identifier/description: Verbose description
* equilibrium/time_slice{:}/boundary/outline/r: Major radius
* equilibrium/time_slice{:}/boundary/outline/z: Height
* equilibrium/time_slice{:}/boundary/lcfs/r: Major radius
* equilibrium/time_slice{:}/boundary/lcfs/z: Height
* equilibrium/time_slice{:}/boundary/geometric_axis/r: Major radius
* equilibrium/time_slice{:}/boundary/geometric_axis/z: Height
* equilibrium/time_slice{:}/boundary_separatrix/outline/r: Major radius
* equilibrium/time_slice{:}/boundary_separatrix/outline/z: Height
* equilibrium/time_slice{:}/boundary_separatrix/geometric_axis/r: Major radius
* equilibrium/time_slice{:}/boundary_separatrix/geometric_axis/z: Height
* equilibrium/time_slice{:}/global_quantities/magnetic_axis/r: Major radius of the magnetic axis
* equilibrium/time_slice{:}/global_quantities/magnetic_axis/z: Height of the magnetic axis
* equilibrium/time_slice{:}/global_quantities/magnetic_axis/b_tor: Total toroidal magnetic field at the magnetic axis
* equilibrium/time_slice{:}/global_quantities/magnetic_axis/b_field_tor: Total toroidal magnetic field at the magnetic axis
* equilibrium/time_slice{:}/global_quantities/magnetic_axis/b_field_phi: Total toroidal magnetic field at the magnetic axis
* equilibrium/time_slice{:}/global_quantities/q_min/value: Minimum q value
* equilibrium/time_slice{:}/global_quantities/q_min/rho_tor_norm: Minimum q position in normalised toroidal flux coordinate
* equilibrium/time_slice{:}/ggd{:}/r{:}/values: One scalar value is provided per element in the grid subset.
* equilibrium/time_slice{:}/ggd{:}/z{:}/values: One scalar value is provided per element in the grid subset.
* equilibrium/time_slice{:}/ggd{:}/psi{:}/values: One scalar value is provided per element in the grid subset.
* equilibrium/time_slice{:}/ggd{:}/phi{:}/values: One scalar value is provided per element in the grid subset.
* equilibrium/time_slice{:}/ggd{:}/theta{:}/values: One scalar value is provided per element in the grid subset.
* equilibrium/time_slice{:}/ggd{:}/j_tor{:}/values: One scalar value is provided per element in the grid subset.
* equilibrium/time_slice{:}/ggd{:}/j_parallel{:}/values: One scalar value is provided per element in the grid subset.
* equilibrium/time_slice{:}/ggd{:}/b_field_r{:}/values: One scalar value is provided per element in the grid subset.
* equilibrium/time_slice{:}/ggd{:}/b_field_tor{:}/values: One scalar value is provided per element in the grid subset.
* equilibrium/time_slice{:}/ggd{:}/b_field_z{:}/values: One scalar value is provided per element in the grid subset.
* equilibrium/time_slice{:}/coordinate_system/grid_type/name: Short string identifier
* equilibrium/time_slice{:}/coordinate_system/grid_type/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* equilibrium/time_slice{:}/coordinate_system/grid_type/description: Verbose description
* equilibrium/time_slice{:}/coordinate_system/grid/dim1: First dimension values
* equilibrium/time_slice{:}/coordinate_system/grid/dim2: Second dimension values
* equilibrium/time_slice{:}/profiles_2d{:}/grid_type/name: Short string identifier
* equilibrium/time_slice{:}/profiles_2d{:}/grid_type/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* equilibrium/time_slice{:}/profiles_2d{:}/grid_type/description: Verbose description
* equilibrium/time_slice{:}/profiles_2d{:}/grid/dim1: First dimension values
* equilibrium/time_slice{:}/profiles_2d{:}/grid/dim2: Second dimension values
* equilibrium/ids_properties/provenance/node{:}/reference{:}/name: Reference name
* equilibrium/ids_properties/provenance/node{:}/reference{:}/timestamp: Date and time (UTC) at which the reference was created, expressed in a human readable form (ISO 8601) : the format of the string shall be : YYYY-MM-DDTHH:MM:SSZ. Example : 2020-07-24T14:19:00Z
* equilibrium/time_slice{:}/ggd{:}/grid/identifier/name: Short string identifier
* equilibrium/time_slice{:}/ggd{:}/grid/identifier/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* equilibrium/time_slice{:}/ggd{:}/grid/identifier/description: Verbose description
* equilibrium/grids_ggd{:}/grid{:}/space{:}/objects_per_dimension{:}/object{:}/geometry: Geometry data associated with the object, its detailed content is defined by ../../geometry_content. Its dimension depends on the type of object, geometry and coordinate considered. 
* equilibrium/time_slice{:}/ggd{:}/grid/space{:}/objects_per_dimension{:}/object{:}/geometry: Geometry data associated with the object, its detailed content is defined by ../../geometry_content. Its dimension depends on the type of object, geometry and coordinate considered. 


IDS plasma_profiles
--------------------

* plasma_profiles/time: Generic time
* plasma_profiles/ids_properties/source: Any comment describing the content of this IDS
* plasma_profiles/ids_properties/provider: User-defined name for this IDS occurrence
* plasma_profiles/ids_properties/creation_date: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* plasma_profiles/profiles_1d{:}/n_i_total_over_n_e: Upper error for "t_i_average"
* plasma_profiles/profiles_1d{:}/momentum_tor: Total ion thermal density (sum over species and charge states)
* plasma_profiles/profiles_1d{:}/zeff: Total plasma toroidal momentum, summed over ion species and electrons weighted by their density and major radius, i.e. sum_over_species(n*R*m*Vphi)
* plasma_profiles/profiles_1d{:}/pressure_ion_total: Upper error for "zeff"
* plasma_profiles/profiles_1d{:}/pressure_thermal: Total (sum over ion species) thermal ion pressure
* plasma_profiles/profiles_1d{:}/pressure_perpendicular: Thermal pressure (electrons+ions)
* plasma_profiles/profiles_1d{:}/pressure_parallel: Total perpendicular pressure (electrons+ions, thermal+non-thermal)
* plasma_profiles/profiles_1d{:}/j_total: Total parallel pressure (electrons+ions, thermal+non-thermal)
* plasma_profiles/profiles_1d{:}/current_parallel_inside: Total parallel current density = average(jtot.B) / B0, where B0 = plasma_profiles/Vacuum_Toroidal_Field/ B0
* plasma_profiles/profiles_1d{:}/j_tor: Parallel current driven inside the flux surface. Cumulative surface integral of j_total
* plasma_profiles/profiles_1d{:}/j_ohmic: Total toroidal current density = average(j_phi/R) / average(1/R)
* plasma_profiles/profiles_1d{:}/j_non_inductive: Ohmic parallel current density = average(J_Ohmic.B) / B0, where B0 = plasma_profiles/Vacuum_Toroidal_Field/ B0
* plasma_profiles/profiles_1d{:}/j_bootstrap: Non-inductive (includes bootstrap) parallel current density = average(jni.B) / B0, where B0 = plasma_profiles/Vacuum_Toroidal_Field/ B0
* plasma_profiles/profiles_1d{:}/conductivity_parallel: Bootstrap current density = average(J_Bootstrap.B) / B0, where B0 = plasma_profiles/Vacuum_Toroidal_Field/ B0
* plasma_profiles/profiles_1d{:}/e_field_parallel: Parallel conductivity
* plasma_profiles/profiles_1d{:}/q: Derivative of the flux surface averaged electrostatic potential with respect to the poloidal flux, multiplied by (-1/2pi). This quantity is the toroidal angular rotation frequency due to the ExB drift, introduced in formula (43) of Hinton and Wong, Physics of Fluids 3082 (1985), also referred to as sonic flow in regimes in which the toroidal velocity is dominant over the poloidal velocity
* plasma_profiles/profiles_1d{:}/magnetic_shear: Safety factor
* plasma_profiles/profiles_1d{:}/time: Magnetic shear, defined as rho_tor/q . dq/drho_tor
* plasma_profiles/global_quantities/ip: Upper error for "pressure_parallel"
* plasma_profiles/global_quantities/current_non_inductive: Total plasma current (toroidal component). Positive sign means anti-clockwise when viewed from above.
* plasma_profiles/global_quantities/current_bootstrap: Total non-inductive current (toroidal component). Positive sign means anti-clockwise when viewed from above.
* plasma_profiles/global_quantities/v_loop: Bootstrap current (toroidal component). Positive sign means anti-clockwise when viewed from above.
* plasma_profiles/global_quantities/li: LCFS loop voltage (positive value drives positive ohmic current that flows anti-clockwise when viewed from above)
* plasma_profiles/global_quantities/li_3: Internal inductance. The li_3 definition is used, i.e. li_3 = 2/R0/mu0^2/Ip^2 * int(Bp^2 dV).
* plasma_profiles/global_quantities/beta_tor: Internal inductance. The li_3 definition is used, i.e. li_3 = 2/R0/mu0^2/Ip^2 * int(Bp^2 dV).
* plasma_profiles/global_quantities/beta_tor_norm: Toroidal beta, defined as the volume-averaged total perpendicular pressure divided by (B0^2/(2*mu0)), i.e. beta_toroidal = 2 mu0 int(p dV) / V / B0^2
* plasma_profiles/global_quantities/beta_pol: Normalised toroidal beta, defined as 100 * beta_tor * a[m] * B0 [T] / ip [MA] 
* plasma_profiles/global_quantities/energy_diamagnetic: Poloidal beta. Defined as betap = 4 int(p dV) / [R_0 * mu_0 * Ip^2]
* plasma_profiles/vacuum_toroidal_field/r0: Time slice of the profiles_1d array used to define the ion composition of the global_quantities/ion array.
* plasma_profiles/vacuum_toroidal_field/b0: Reference major radius where the vacuum toroidal magnetic field is given (usually a fixed position such as the middle of the vessel at the equatorial midplane)
* plasma_profiles/code/name: Upper error for "time_width"
* plasma_profiles/code/commit: Index in the error_description list for "time_width"
* plasma_profiles/code/version: Time
* plasma_profiles/code/repository: Name of software generating IDS
* plasma_profiles/code/parameters: Short description of the software (type, purpose)
* plasma_profiles/code/output_flag: Unique commit reference of software


IDS plasma_transport
---------------------

* plasma_transport/time: Generic time
* plasma_transport/ids_properties/homogeneous_time: Any comment describing the content of this IDS
* plasma_transport/ids_properties/source: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* plasma_transport/ids_properties/creation_date: Source of the data (any comment describing the origin of the data : code, path to diagnostic signals, processing method, ...). Superseeded by the new provenance structure.
* plasma_transport/vacuum_toroidal_field/b0: Lower error for "r0"
* plasma_transport/model{:}/flux_multiplier: Any comment describing the model
* plasma_transport/code/repository: Name of software generating IDS
* plasma_transport/code/parameters: Short description of the software (type, purpose)
* plasma_transport/code/output_flag: Unique commit reference of software


IDS plasma_sources
-------------------

* plasma_sources/time: Generic time
* plasma_sources/ids_properties/homogeneous_time: Any comment describing the content of this IDS
* plasma_sources/ids_properties/source: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* plasma_sources/ids_properties/creation_date: Source of the data (any comment describing the origin of the data : code, path to diagnostic signals, processing method, ...). Superseeded by the new provenance structure.
* plasma_sources/vacuum_toroidal_field/b0: Lower error for "r0"
* plasma_sources/code/commit: Name of software generating IDS
* plasma_sources/code/version: Short description of the software (type, purpose)
* plasma_sources/code/repository: Unique commit reference of software
* plasma_sources/code/parameters: Unique version (tag) of software
* plasma_sources/code/output_flag: URL of software repository


IDS edge_profiles
------------------

* edge_profiles/time: Generic time
* edge_profiles/ids_properties/comment: Any comment describing the content of this IDS
* edge_profiles/ids_properties/homogeneous_time: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* edge_profiles/grid_ggd{:}/time: Time
* edge_profiles/ggd{:}/time: Time
* edge_profiles/code/name: Name of software generating IDS
* edge_profiles/code/commit: Unique commit reference of software
* edge_profiles/code/version: Unique version (tag) of software
* edge_profiles/code/repository: URL of software repository
* edge_profiles/code/parameters: List of the code specific parameters in XML format
* edge_profiles/code/output_flag: Output flag : 0 means the run is successful, other values mean some difficulty has been encountered, the exact meaning is then code specific. Negative values mean the result shall not be used.
* edge_profiles/ggd{:}/zeff{:}/grid_index: Index of the grid used to represent this quantity
* edge_profiles/ggd{:}/zeff{:}/grid_subset_index: Index of the grid subset the data is provided on. Corresponds to the index used in the grid subset definition: grid_subset(:)/identifier/index
* edge_profiles/ggd{:}/zeff{:}/values: One scalar value is provided per element in the grid subset.
* edge_profiles/ggd{:}/electrons/temperature{:}/grid_index: Index of the grid used to represent this quantity
* edge_profiles/ggd{:}/electrons/temperature{:}/grid_subset_index: Index of the grid subset the data is provided on. Corresponds to the index used in the grid subset definition: grid_subset(:)/identifier/index
* edge_profiles/ggd{:}/electrons/temperature{:}/values: One scalar value is provided per element in the grid subset.
* edge_profiles/ggd{:}/electrons/density{:}/grid_index: Index of the grid used to represent this quantity
* edge_profiles/ggd{:}/electrons/density{:}/grid_subset_index: Index of the grid subset the data is provided on. Corresponds to the index used in the grid subset definition: grid_subset(:)/identifier/index
* edge_profiles/ggd{:}/electrons/density{:}/values: One scalar value is provided per element in the grid subset.


IDS edge_transport
-------------------

* edge_transport/time: Generic time
* edge_transport/ids_properties/comment: Any comment describing the content of this IDS
* edge_transport/ids_properties/homogeneous_time: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* edge_transport/grid_ggd{:}/time: Time
* edge_transport/code/name: Name of software generating IDS
* edge_transport/code/commit: Unique commit reference of software
* edge_transport/code/version: Unique version (tag) of software
* edge_transport/code/repository: URL of software repository
* edge_transport/code/parameters: List of the code specific parameters in XML format
* edge_transport/code/output_flag: Output flag : 0 means the run is successful, other values mean some difficulty has been encountered, the exact meaning is then code specific. Negative values mean the result shall not be used.
* edge_transport/model{:}/identifier/name: Short string identifier
* edge_transport/model{:}/identifier/index: Integer identifier (enumeration index within a list). Private identifier values must be indicated by a negative index.
* edge_transport/model{:}/identifier/description: Verbose description
* edge_transport/model{:}/ggd{:}/time: Time
* edge_transport/model{:}/ggd_fast{:}/power{:}/grid_index: Index of the grid used to represent this quantity
* edge_transport/model{:}/ggd_fast{:}/power{:}/grid_subset_index: Index of the grid subset the data is provided on. Corresponds to the index used in the grid subset definition: grid_subset(:)/identifier/index
* edge_transport/model{:}/ggd_fast{:}/power{:}/value: Scalar value of the quantity on the grid subset (corresponding to a single local position or to an integrated value over the subset)
* edge_transport/model{:}/ggd{:}/electrons/particles/flux{:}/grid_index: Index of the grid used to represent this quantity
* edge_transport/model{:}/ggd{:}/electrons/particles/flux{:}/grid_subset_index: Index of the grid subset the data is provided on. Corresponds to the index used in the grid subset definition: grid_subset(:)/identifier/index
* edge_transport/model{:}/ggd{:}/electrons/particles/flux{:}/values: One scalar value is provided per element in the grid subset.
* edge_transport/model{:}/ggd{:}/electrons/energy/flux{:}/grid_index: Index of the grid used to represent this quantity
* edge_transport/model{:}/ggd{:}/electrons/energy/flux{:}/grid_subset_index: Index of the grid subset the data is provided on. Corresponds to the index used in the grid subset definition: grid_subset(:)/identifier/index
* edge_transport/model{:}/ggd{:}/electrons/energy/flux{:}/values: One scalar value is provided per element in the grid subset.


IDS plasma_profiles_1
----------------------

* plasma_profiles_1/time: Generic time
* plasma_profiles_1/ids_properties/source: Any comment describing the content of this IDS
* plasma_profiles_1/ids_properties/provider: User-defined name for this IDS occurrence
* plasma_profiles_1/ids_properties/creation_date: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* plasma_profiles_1/profiles_1d{:}/n_i_total_over_n_e: Upper error for "t_i_average"
* plasma_profiles_1/profiles_1d{:}/momentum_tor: Total ion thermal density (sum over species and charge states)
* plasma_profiles_1/profiles_1d{:}/zeff: Total plasma toroidal momentum, summed over ion species and electrons weighted by their density and major radius, i.e. sum_over_species(n*R*m*Vphi)
* plasma_profiles_1/profiles_1d{:}/pressure_ion_total: Upper error for "zeff"
* plasma_profiles_1/profiles_1d{:}/pressure_thermal: Total (sum over ion species) thermal ion pressure
* plasma_profiles_1/profiles_1d{:}/pressure_perpendicular: Thermal pressure (electrons+ions)
* plasma_profiles_1/profiles_1d{:}/pressure_parallel: Total perpendicular pressure (electrons+ions, thermal+non-thermal)
* plasma_profiles_1/profiles_1d{:}/j_total: Total parallel pressure (electrons+ions, thermal+non-thermal)
* plasma_profiles_1/profiles_1d{:}/current_parallel_inside: Total parallel current density = average(jtot.B) / B0, where B0 = plasma_profiles/Vacuum_Toroidal_Field/ B0
* plasma_profiles_1/profiles_1d{:}/j_tor: Parallel current driven inside the flux surface. Cumulative surface integral of j_total
* plasma_profiles_1/profiles_1d{:}/j_ohmic: Total toroidal current density = average(j_phi/R) / average(1/R)
* plasma_profiles_1/profiles_1d{:}/j_non_inductive: Ohmic parallel current density = average(J_Ohmic.B) / B0, where B0 = plasma_profiles/Vacuum_Toroidal_Field/ B0
* plasma_profiles_1/profiles_1d{:}/j_bootstrap: Non-inductive (includes bootstrap) parallel current density = average(jni.B) / B0, where B0 = plasma_profiles/Vacuum_Toroidal_Field/ B0
* plasma_profiles_1/profiles_1d{:}/conductivity_parallel: Bootstrap current density = average(J_Bootstrap.B) / B0, where B0 = plasma_profiles/Vacuum_Toroidal_Field/ B0
* plasma_profiles_1/profiles_1d{:}/e_field_parallel: Parallel conductivity
* plasma_profiles_1/profiles_1d{:}/q: Derivative of the flux surface averaged electrostatic potential with respect to the poloidal flux, multiplied by (-1/2pi). This quantity is the toroidal angular rotation frequency due to the ExB drift, introduced in formula (43) of Hinton and Wong, Physics of Fluids 3082 (1985), also referred to as sonic flow in regimes in which the toroidal velocity is dominant over the poloidal velocity
* plasma_profiles_1/profiles_1d{:}/magnetic_shear: Safety factor
* plasma_profiles_1/profiles_1d{:}/time: Magnetic shear, defined as rho_tor/q . dq/drho_tor
* plasma_profiles_1/global_quantities/ip: Upper error for "pressure_parallel"
* plasma_profiles_1/global_quantities/current_non_inductive: Total plasma current (toroidal component). Positive sign means anti-clockwise when viewed from above.
* plasma_profiles_1/global_quantities/current_bootstrap: Total non-inductive current (toroidal component). Positive sign means anti-clockwise when viewed from above.
* plasma_profiles_1/global_quantities/v_loop: Bootstrap current (toroidal component). Positive sign means anti-clockwise when viewed from above.
* plasma_profiles_1/global_quantities/li: LCFS loop voltage (positive value drives positive ohmic current that flows anti-clockwise when viewed from above)
* plasma_profiles_1/global_quantities/li_3: Internal inductance. The li_3 definition is used, i.e. li_3 = 2/R0/mu0^2/Ip^2 * int(Bp^2 dV).
* plasma_profiles_1/global_quantities/beta_tor: Internal inductance. The li_3 definition is used, i.e. li_3 = 2/R0/mu0^2/Ip^2 * int(Bp^2 dV).
* plasma_profiles_1/global_quantities/beta_tor_norm: Toroidal beta, defined as the volume-averaged total perpendicular pressure divided by (B0^2/(2*mu0)), i.e. beta_toroidal = 2 mu0 int(p dV) / V / B0^2
* plasma_profiles_1/global_quantities/beta_pol: Normalised toroidal beta, defined as 100 * beta_tor * a[m] * B0 [T] / ip [MA] 
* plasma_profiles_1/global_quantities/energy_diamagnetic: Poloidal beta. Defined as betap = 4 int(p dV) / [R_0 * mu_0 * Ip^2]
* plasma_profiles_1/vacuum_toroidal_field/r0: Time slice of the profiles_1d array used to define the ion composition of the global_quantities/ion array.
* plasma_profiles_1/vacuum_toroidal_field/b0: Reference major radius where the vacuum toroidal magnetic field is given (usually a fixed position such as the middle of the vessel at the equatorial midplane)
* plasma_profiles_1/code/name: Upper error for "time_width"
* plasma_profiles_1/code/commit: Index in the error_description list for "time_width"
* plasma_profiles_1/code/version: Time
* plasma_profiles_1/code/repository: Name of software generating IDS
* plasma_profiles_1/code/parameters: Short description of the software (type, purpose)
* plasma_profiles_1/code/output_flag: Unique commit reference of software


IDS plasma_transport_1
-----------------------

* plasma_transport_1/time: Generic time
* plasma_transport_1/ids_properties/homogeneous_time: Any comment describing the content of this IDS
* plasma_transport_1/ids_properties/source: This node must be filled (with 0, 1, or 2) for the IDS to be valid. If 1, the time of this IDS is homogeneous, i.e. the time values for this IDS are stored in the time node just below the root of this IDS. If 0, the time values are stored in the various time fields at lower levels in the tree. In the case only constant or static nodes are filled within the IDS, homogeneous_time must be set to 2
* plasma_transport_1/ids_properties/creation_date: Source of the data (any comment describing the origin of the data : code, path to diagnostic signals, processing method, ...). Superseeded by the new provenance structure.
* plasma_transport_1/vacuum_toroidal_field/b0: Lower error for "r0"
* plasma_transport_1/model{:}/flux_multiplier: Any comment describing the model
* plasma_transport_1/code/repository: Name of software generating IDS
* plasma_transport_1/code/parameters: Short description of the software (type, purpose)
* plasma_transport_1/code/output_flag: Unique commit reference of software


