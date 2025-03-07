#!/Applications/Anaconda/python37/bin/python3.7
# coding: utf8

# Built-in
import os
import argparse
import warnings
import traceback


import scipy.io as scpio


# Common
# tofu plugin
import tofu_west as tfw


_diag = ['interf', 'bolo']
_shot = 0
_run = 0
_occ = 0
_usr = 'imas_public'
_machine = 'west'
_save_mode = 'ids'
_save_path = './'
_save_name = None


def call_tfw_synthdiag(diag=_diag, shot=_shot, run=_run, occ=_occ,
                       usr=_usr, machine=_machine, input_file=None,
                       shot_write=None, run_write=None, occ_write=None,
                       usr_write=None, machine_write=None,
                       save_mode=_save_mode, save_path=_save_path,
                       save_name=_save_name):
    msg = "You must specify a valid diag (or list of such) !\n"
    msg += "Not valid: "+str(diag)
    assert type(diag) in [str,list] and not diag=='' and not diag==[], msg
    msg = "Args usr, machine must be str !"
    assert all([type(ss) is str for ss in [usr,machine]]), msg
    msg = "Args shot, run and occ must be int !"
    assert all([type(ss) is int for ss in [shot,run, occ]]), msg
    assert save_mode in ['ids','mat']
    assert type(save_path) is str

    # Preformat inputs and def param
    if type(diag) is str:
        diag = [diag]
    for ii in range(0,len(diag)):
        diag[ii] = diag[ii].lower()

    dfunc = {'sxr':tfw.SXR, 'bolo':tfw.Bolo,
             'interf':tfw.Interfero, 'brem':tfw.SpectroVisZeff}
    dtofu = {'sxr':'SXR', 'bolo':'Bolo', 'interf':'Interfero',
             'brem':'SpectroVisZeff'}

    # Set imas parameters for reading input data
    dimas_read = {'shot':shot, 'run':run, 'occ':occ,
                  'usr':usr, 'machine':machine}

    if save_mode=='ids':
        # Set imas parameters for writing output signal
        shot_write = shot if shot_write is None else shot_write
        run_write = run if run_write is None else run_write
        occ_write = occ if occ_write is None else occ_write
        usr_write = usr if usr_write is None else usr_write
        machine_write = machine if machine_write is None else usr_write
        dimas_write = {'shot':shot_write, 'run':run_write, 'occ':occ_write,
                       'usr':usr_write, 'machine':machine_write}

    # Launch diags
    msg0 = "Could not launch "
    for ii in range(0,len(diag)):
        if diag[ii] in dfunc.keys():
            try:
                print("")
                print("Trying synthetic {0}".format(dtofu[diag[ii]]))
                out, kh = dfunc[diag[ii]].calc_signal(dimas_read=dimas_read,
                                                      input_file=input_file,
                                                      dimas_write=None,
                                                      plot=False)

                if save_mode=='ids':
                    try:
                        out.save_to_imas(**dimas_write)
                    except Exception as err:
                        msg = str(err)
                        msg += "\nCould not save computed synthetic signal to:\n"
                        msg += "out.save_to_imas(**dimas_write), with:\n"
                        msg += "dimas_write = {"
                        msg += ", ".join([kk+":{0}".format(dimas_write[kk])
                                          for kk in dimas_write.keys()])+"}\n"
                        msg += "\nInput parameters were:\n"
                        msg += "dimas_read = "
                        msg += ", ".join([kk+":{0}".format(dimas_write[kk])
                                          for kk in dimas_write.keys()])+"}\n"
                        msg += "input_file = %s"%input_file
                        warnings.warn(msg)
                else:
                    pfe = None
                    try:
                        # Format output dictionnary to be saved 
                        dout = {
                            'shot': shot,
                            't': out.t,
                            'data': out.data,
                            'units_t': out.dlabels['t']['units'],
                            'units_data': out.dlabels['data']['units'],
                            'channels': out.dchans('Name'),
                            'tofu_west_version': tfw.__version__,
                        }

                        # Save to specified path + filename + extension
                        if save_name is None:
                            save_name = 'tofu_{0}_{1}.mat'.format(shot,diag[ii])
                        assert type(save_name) is str
                        if save_name[-4:]!='.mat':
                            assert len(save_name.split('.'))==1
                            save_name += '.mat'
                        pfe = os.path.join(os.path.abspath(save_path),
                                           save_name)
                        scpio.savemat(pfe, dout)
                    except Exception as err:
                        msg = str(traceback.format_exc())
                        msg += "\n\n" + str(err)
                        msg += "\nCould not save computed synthetic signal to:\n"
                        msg += "scpio.savemat({0}, dout)".format(pfe)
                        warnings.warn(msg)


            except Exception as err:
                msg = str(traceback.format_exc())
                msg += "\n\n" + str(err)
                msg += "\nCould not compute signal:"
                msg += "tfw.{0}.calc_signal(dimas_read), with:\n".format(dtofu[diag[ii]])
                msg += "dimas_read = {"
                msg += ", ".join([kk+":{0}".format(dimas_read[kk])
                                  for kk in dimas_read.keys()]) + "}"
                warnings.warn(msg)




if __name__ == '__main__':

    # Parse input arguments
        msg = 'Wrapper to launch synthetic diags of tofu_west from bash'
        parser = argparse.ArgumentParser(description = msg)

        msg = "diag, can be:"
        msg += "\n    - 'bolo':         Bolometry"
        msg += "\n    - 'sxr':          Soft X-Rays"
        msg += "\n    - 'interf':       Interferometry"
        msg += "\n    - 'brem':         Bremsshtrahlung visible spectroscopy"
        parser.add_argument('-diag', type=str, help=msg, default=_diag)

        parser.add_argument('-shot', type=int, help='shot', default=_shot)
        parser.add_argument('-run', type=int, help='run', default=_run)
        parser.add_argument('-occ', type=int, help='occ', default=_occ)
        parser.add_argument('-usr', type=str, help='user', default=_usr)
        parser.add_argument('-machine', type=str, help='machine',
                            default=_machine)
        parser.add_argument('-input_file', type=str, help='.mat for brem',
                            default=None)
        parser.add_argument('-save_mode', type=str, help='ids or mat',
                            default=_save_mode)
        parser.add_argument('-save_path', type=str, help='for mat only',
                            default=_save_path)
        parser.add_argument('-save_name', type=str, help='for mat only',
                            default=_save_name)

        args = parser.parse_args()

        # Call wrapper function
        call_tfw_synthdiag(diag=args.diag, shot=args.shot, run=args.run,
                           occ=args.occ, usr=args.usr, machine=args.machine,
                           save_mode=args.save_mode, save_path=args.save_path,
                           save_name=args.save_name, input_file=args.input_file)


