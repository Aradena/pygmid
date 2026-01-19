import subprocess
import logging
from paramiko import SSHClient, SSHException


class SpectreSimulator:
    def __init__(self, _config, *args):
        self.__args = list(args)
        self._config = _config
        self._init_ssh()

    @property
    def directory(self):
        return self.__args[-1]
    
    @directory.setter
    def directory(self, dir):
        self.__args[-1] = dir

    def run(self, filename: str):
        infile = filename
        try:
            cmd_args = ['spectre', filename] + [*self.__args]
            cp = subprocess.run(cmd_args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            logging.info(f"Error executing process\n\n{e}")
            return
        
        return self.__args[-1]

    def run_through_ssh(self, filename: str):

        cmd_args = self._config['SSH']['LOAD_PDK']+" && "+"cd "+self._config['SSH']['RUNDIR'] +' && '+ str(['spectre', filename] + [*self.__args])
        try:
            stdin, stdout, stderr = self.client.exec_command(cmd_args)
            print(stdout.readlines(), stderr.readlines())

        except SSHException as e:
            logging.info(f"Error executing command\n\n{e}")
            return

        return self.__args[-1]

    def _init_ssh(self):
        if 'SSH' in self._config.keys():
            self.client = SSHClient()
            self.client.load_system_host_keys()
            self.client.connect(hostname=self._config['SSH']['HOSTNAME'],
                                username=self._config['SSH']['USERNAME'],
                                password=self._config['SSH']['KEY_FILENAME'])

            # stdin, stdout, stderr = self.client.exec_command('pwd')
            # print(stdout.readlines())
            # each time exec_command is executed, back to ~/


