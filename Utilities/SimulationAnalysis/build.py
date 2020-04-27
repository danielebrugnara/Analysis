#!/usr/bin/env python3
import os
import argparse

def clean():
    os.system("rm -f *.pcm")
    os.system("rm -f Makefile")
    os.system("rm -f *.a")
    os.system("rm -f *.rootmap")
    os.system("rm -f G__*")
    os.system("rm -f cmake_install.cmake")
    os.system("rm -f -r CMakeFiles")
    os.system("rm -f CMakeCache.txt")
    os.system("rm -f main")

def build():
    os.system("cmake .")
    os.system("make")

def main():
    parser = argparse.ArgumentParser(description='Run and build analysisn simulation')
    parser.add_argument('--clean',
                    action='store_true',
                    dest='clean',
                    help='cleans previous build')
    parser.add_argument('--build',
                    action='store_true',
                    dest='build',
                    help='builds')
    parser.add_argument('--clean_and_build',
                    action='store_true',
                    dest='clean_and_build',
                    help='cleans and builds previous build')

    args = parser.parse_args()

    if (args.clean):
        clean()
        return 0
    if (args.build):
        build()
        return 0
    if (args.clean_and_build):
        build()
        clean()
        return 0

    


if __name__ == '__main__':
    main()
