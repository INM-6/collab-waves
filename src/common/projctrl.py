'''
A collection of routines to facilitate handling a project that produces data
and figures.
'''
import inspect
import os.path
import shutil

import matplotlib.pyplot as plt

__updated__ = "2016-02-15"

if __name__ == '__main__':
    pass


class ProjectControl(object):
    '''
    This class manages paths to results and figures created for a specific
    task in a project.

    Here, a project task refers to the following:
        * an (analysis) script in python
        * a project directory, or workspace, containing all related scripts to
          the project
        * a project data directory containing all results and figures created
          in the project

    By default, this class determines the name of the current running script
    x.py. Next, the script creates, if necessary, corresponding results and
    figures directories in results/x and figs/x located within the project data
    directory.

    This class gives access to the paths to results and figures, and figures
    can be directly saved to a variety of formats in the correct location.

    The constructor will determine project and script name, and create (if
    required) directories to store results and figures.

    Parameters:
    -----------
    project_name : string
        Name of the project. If set to None, the project name is
        guessed from the directory the calling script is in coming from
        workspace_base.
        Default: None.
    script_name :string
        Name of the script. If set to None, the script name is guessed
        from the name of the calling script.
        Default: None.
    workspace_base : string
        This is the base directory of the workspace. Projects are
        assumed to be subdirectories in this directory.
        Default: 'Projects'
    data_base : string
        This is the base directory of result and figure files. Projects are
        assumed to be subdirectories in this directory.
        Default: 'ProjectsData'
    clear_data : boolean
        If True, both result and figure directories are emptied (i.e.,
        all files are deleted). Use with care!
        Default: False

    Attributes:
    -----------
    result_path : string
        Path where to store results, e.g., a file might be located in
        data_base/results/x/01.dat where x is script_name.
    figure_path : string
        Top level path where figures are stored. Actual figures are in
        subdirectories corresponding to file extensions, e.g., a file might
        be located in data_base/figs/x/jpg/01.jpg where x is script_name.
    result_domain : string
        Parent directory where all project results are stored, e.g.,
        data_base/results/.
    figure_domain : string
        Parent directory where all project figures are stored, e.g.,
        data_base/figs/.
    '''

    def __init__(
            self, project_name=None, script_name=None,
            workspace_base='Projects', data_base='ProjectsData',
            clear_data=False):
        '''
        Constructor.
        '''

        # get name of current script and its path -- this is the file calling
        # this function
        (callingdir, callingfile) = os.path.split(inspect.stack()[-1][1])

        # if no project name is given, extract it as the folder following the
        # folder given by the variable workspace_base. Example: in
        # '/user/foobar/Project/myproject' the project name becomes 'myproject'
        # when using the default workspace_base 'Project'.
        if project_name is None:
            # disect the calling directory into its individual parts of the
            # path
            folders = []
            path = callingdir
            while 1:
                path, folder = os.path.split(path)

                if folder != workspace_base and folder != '':
                    folders.append(folder)
                else:
                    if folder == '':
                        raise ValueError('Error: could not detect project \
                            name during initialization')

                    break

            # select last folder added as project name
            project_name = folders[-1]

        # if no script name is specified, extract it from the name of the
        # calling python script
        if script_name is None:
            script_name = os.path.splitext(callingfile)[0]
            if script_name == '':
                raise ValueError('Error: could not detect script name \
                    during initialization')

        # build complete paths of figures and results in the user's home
        # directory
        self.result_domain = os.path.expanduser(
            '~' + os.path.sep + data_base +
            os.path.sep + project_name + os.path.sep + 'results')
        self.result_path = os.path.expanduser(
            self.result_domain + os.path.sep + script_name + os.path.sep)
        self.figure_domain = os.path.expanduser(
            '~' + os.path.sep + data_base +
            os.path.sep + project_name + os.path.sep + 'figs')
        self.figure_path = os.path.expanduser(
            self.figure_domain + os.path.sep + script_name + os.path.sep)

        # create paths if not yet existent, empty them if clear_data==True
        if clear_data:
            if os.path.exists(self.result_path):
                shutil.rmtree(self.result_path)
            if os.path.exists(self.figure_path):
                shutil.rmtree(self.figure_path)

        if not os.path.exists(self.result_path):
            try:
                os.makedirs(self.result_path)
            except OSError:
                pass

        if not os.path.exists(self.figure_path):
            os.makedirs(self.figure_path + os.path.sep + 'jpg')
            os.makedirs(self.figure_path + os.path.sep + 'png')
            os.makedirs(self.figure_path + os.path.sep + 'eps')
            os.makedirs(self.figure_path + os.path.sep + 'pdf')
            os.makedirs(self.figure_path + os.path.sep + 'mov')

    def get_paths(self):
        '''
        Returns a tuple of two elements, the first containing the base path for
        results, the second the base path for figures of the project.

        Returns:
        --------
        Tuple
            Tuple (r,f), where r is a string containing the results path
            and f is a string containing the figure path.
        '''
        return (self.result_path, self.figure_path)

    def save_fig(
            self, filename,
            save_png=True, save_jpg=False, save_eps=True, save_pdf=False,
            dpi=150):
        '''
        Save the current figure to files in the project's figure directory.

        Parameters:
        -----------
        filename (string):
            Base filename of the files (without path and without
            extension).
        save_png (bool):
            If True, a PNG version of the figure is saved. Default: True
        save_jpg (bool):
            If True, a JPG version of the figure is saved. Default: False
        save_eps (bool):
            If True, a EPS version of the figure is saved. Default: True
        save_pdf (bool):
            If True, a PDF version of the figure is saved. Default: False
        dpi (int):
            DPI resolution. Default: 150.
        '''

        # strip any path information
        filename = os.path.split(filename)[1]

        # save figures
        if save_jpg:
            plt.savefig(
                self.figure_path + 'jpg' + os.path.sep + filename + '.jpg',
                dpi=dpi)
        if save_png:
            plt.savefig(
                self.figure_path + 'png' + os.path.sep + filename + '.png',
                dpi=dpi)
        if save_eps:
            plt.savefig(
                self.figure_path + 'eps' + os.path.sep + filename + '.eps',
                dpi=dpi)
        if save_pdf:
            plt.savefig(
                self.figure_path + 'pdf' + os.path.sep + filename + '.pdf',
                dpi=dpi)
