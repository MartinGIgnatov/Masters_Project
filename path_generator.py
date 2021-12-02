# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 12:46:52 2021

@author: marti
"""

import os

def generate_path(rel_path, file_name, extension):
    """
    generate file path depending on OS
    rel_path   ---   defines the relative path, only with folders. can be a str a list of str or None for no folders.
    file_name  ---   name of file
    extension     ---   extension of file
    """
    if not isinstance(extension, str):
        raise Exception("Wrong data type for extension")
    if not isinstance(file_name, str):
        raise Exception("Wrong data type for file_name")
    if not isinstance(rel_path, (str,list,None)):
        raise Exception("Wrong data type for rel_path")
    if isinstance(rel_path, list):
        for el in rel_path:
            if not isinstance(el, str):
                raise Exception("Wrong data type for rel_path element")
    
    if os.name == "nt":
        path = ""
        if isinstance(rel_path, str):
            path =  f".\\{rel_path}\\{file_name}.{extension}"
        elif isinstance(rel_path, None):
            path =  f".\\{file_name}.{extension}"
        else:
            path = ".\\"
            for el in rel_path:
                path = path + el + "\\"
            path = path + f"{file_name}.{extension}"
        return path
            
    elif os.name == "posix":   #### This is wrong
        if isinstance(rel_path, str):
            path =  f"{rel_path}/{file_name}.{extension}"
        elif isinstance(rel_path, None):
            path =  f"{file_name}.{extension}"
        else:
            path = ""
            for el in rel_path:
                path = path + el + ""
            path = path + f"{file_name}.{extension}"
        return path
    else:
        raise Exception("This OS is not recognized.")