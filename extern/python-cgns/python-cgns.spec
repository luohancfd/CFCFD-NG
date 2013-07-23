# To activate HDF5 support use "--with hdf5"
# Support for HDF5 is deactivated by default
%bcond_with hdf5

# Set Low Level Library name
%if %{with hdf5}
%define lllib hdf5
%else
%define lllib adf
%endif

# Disable debuginfo
%define debug_package %{nil}

# Python modules install locations
%define python_sitelib  %(%{__python} -c "from distutils.sysconfig import get_python_lib; print get_python_lib()")
%define python_sitearch %(%{__python} -c "from distutils.sysconfig import get_python_lib; print get_python_lib(1)")

# CGNS installed version
%define CGNSversion %(rpm -q --qf '%{VERSION}' cgns-%{lllib})

Summary: Python bindings for CGNS library
Name: python-cgns
Version: 0.3
Release: 2%{?dist}.fluorem
License: GPLv3
Vendor:  Fluorem SAS
Group: Development/Languages
URL: http://pypi.python.org/pypi/CGNS
Source0: http://pypi.python.org/packages/source/C/CGNS/%{name}-%{version}.tar.gz
BuildRequires: cgns-%{lllib},cgns-%{lllib}-devel
BuildRequires: python,python-devel
BuildRequires: swig
%if %{with hdf5}
BuildRequires: hdf5 => 1.8, hdf5-devel => 1.8
%endif
BuildRoot: %{_tmppath}/%{name}-%{lllib}-%{version}-%{release}-root

%description
Python CGNS is the Python language binding for the CGNS library.
It consists of several Python Extensions, which provides methods for
reading, writing and manipulation of CGNS files.

%package -n %{name}%{CGNSversion}-%{lllib}
Summary: Python bindings for CGNS library
License: GPLv3
Group: Development/Languages
Requires: python-abi = %(%{__python} -c "import sys ; print sys.version[:3]")
Requires: libcgns = %{CGNSversion}
%if %{with hdf5}
Requires: hdf5 => 1.8
Conflicts: %{name}-adf
Conflicts: cgns-adf
%else
Conflicts: %{name}-hdf5
Conflicts: cgns-hdf5
%endif

%description -n %{name}%{CGNSversion}-%{lllib}
Python CGNS is the Python language binding for the CGNS library.
It consists of several Python Extensions, which provides methods for
reading, writing and manipulation of CGNS files.

%prep
%setup -q

%build
# Copy CGNS lib include file to the build directory
cp -a %{_libdir}/cgns/include/cgnslib.h src/
if [ -f %{_libdir}/cgns/include/cgnstypes.h ]; then
  cp -a %{_libdir}/cgns/include/cgnstypes.h src/
fi

# Change the name of "const char *" arguments in CGNS lib include file
# This is to ensure that cstring_bounded_output SWIG macros will not
# affect input only arguments that are declared as "const char *".
cp -a src/cgnslib.h src/cgnslib.h.fixconstchar
sed -e '/char[[:space:]][[:space:]]*const[[:space:]][[:space:]]*\*[[:space:]]*[A-Za-z0-9_][A-Za-z0-9_]*[,)]/ s/char[[:space:]][[:space:]]*const[[:space:]][[:space:]]*\*[[:space:]]*/char const *varconst_/g' \
    -e '/const[[:space:]][[:space:]]*char[[:space:]][[:space:]]*\*[[:space:]]*[A-Za-z0-9_][A-Za-z0-9_]*[,)]/ s/const[[:space:]][[:space:]]*char[[:space:]][[:space:]]*\*[[:space:]]*/const char *varconst_/g' src/cgnslib.h.fixconstchar > src/cgnslib.h

# Fix to paths in setup.py config file
cp -a setup.py setup.py.fixpath
sed -e 's@swig_opts=\["-I/usr/include"\]@swig_opts=["-I/usr/include -v"], include_dirs=["/usr/include"], library_dirs=["'%{_libdir}'/cgns/lib"]@' setup.py.fixpath > setup.py

# Build
env CFLAGS="$RPM_OPT_FLAGS" %{__python} setup.py build

%install
%{__rm} -rf $RPM_BUILD_ROOT
%{__python} setup.py install --root=$RPM_BUILD_ROOT

%clean
%{__rm} -rf $RPM_BUILD_ROOT

%files -n %{name}%{CGNSversion}-%{lllib}
%defattr(-,root,root,-)
%dir %{python_sitearch}/CGNS
%dir %{_datadir}/doc/%{name}-%{version}
%dir %{_datadir}/doc/%{name}-%{version}/examples
%{python_sitearch}/CGNS/*.py*
%{python_sitearch}/CGNS/*.so
%{_datadir}/doc/%{name}-%{version}/LICENSE
%{_datadir}/doc/%{name}-%{version}/examples/*

%changelog
* Thu Mar 28 2011 Lionel Gamet <Lionel.Gamet@fluorem.com>
- Updated to support CGNS 3.1
- Added CGNS version number in binary package name
* Thu Mar 03 2011 Lionel Gamet <Lionel.Gamet@fluorem.com>
- Updated to 0.3
* Thu Jan 07 2011 Lionel Gamet <Lionel.Gamet@fluorem.com>
- Corrections on const char * handling
- Revisited examples
* Thu Jan 06 2011 Lionel Gamet <Lionel.Gamet@fluorem.com>
- Added optional command line support for HDF5
- Distinct naming convention if compiling against cgns-hdf5 or cgns-adf
* Thu Dec 23 2010 Lionel Gamet <Lionel.Gamet@fluorem.com>
- Changed versioning to 0.2
* Thu Nov 24 2010 Lionel Gamet <Lionel.Gamet@fluorem.com>
- Initial build.
