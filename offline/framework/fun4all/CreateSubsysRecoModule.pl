#!/usr/bin/perl

use strict;
use Getopt::Long;

sub CreateInclude;
sub CreateImplementation;
sub CreateMakefile;
sub CreateAutogen;
sub CreateConfigure;

if ($#ARGV < 0)
{
    print "Usage:\n";
    print "CreateSubsysRecoModule.pl <Module Name>\n";
    print "options:\n";
    print "--all : create also autogen.sh, configure.ac and Makefile.am\n";
    print "--overwrite : overwrite existing files\n";
    exit(0);
}

my $createall;
my $overwrite;
GetOptions("all" => \$createall, "overwrite" => \$overwrite);
my $modulename = $ARGV[0];
my %listoffiles = ();
my $classname = sprintf("%s",$ARGV[0]);
if ($classname =~ m/[^a-zA-Z0-9_]/)
{
    print "Module name contains invalid characters - allowed are alphanumeric (lower and upper caps) and _\n";
    exit(1);
}
if (substr($classname,0,1) =~ m/[0-9_]/)
{
    print "Module name must start with a letter\n";
    exit(1);
}

my $includefile = sprintf("%s.h", $classname);
my $implementfile = sprintf("%s.cc", $classname);
$listoffiles{$includefile} = $classname;
$listoffiles{$implementfile} = $classname;
if (defined $createall)
{
    $listoffiles{"autogen.sh"} = $classname;
    $listoffiles{"configure.ac"} = $classname;
    $listoffiles{"Makefile.am"} = $classname;
}

# check if files exist if overwrite is not set

if (! defined $overwrite)
{
    foreach my $file (keys %listoffiles)
    {
	if (-f $file)
	{
	    print "$file exists but overwrite option not set\n";
	    exit(1);
	}

    }
}

CreateInclude($includefile);
CreateImplementation($implementfile);
if (defined $createall)
{
    CreateAutogen();
    CreateMakefile();
     CreateConfigure();
}
exit(0);

sub CreateImplementation()
{
    my $file = shift;
    open(F,">$file");
    print F "#include \"$classname.h\"\n";
    print F "\n";

    print F "#include <fun4all/Fun4AllReturnCodes.h>\n";
    print F "#include <phool/PHCompositeNode.h>\n";
    print F "\n";

    print F "using namespace std;\n";
    print F "\n";

    print F "$classname\:\:$classname(const std::string &name):\n";
    print F " SubsysReco(name)\n";
    print F "{\n";
    print F "  cout << \"Calling ctor\" << endl;\n";
    print F "}\n";
    print F "\n";

    print F "$classname\:\:~$classname()\n";
    print F "{\n";
    print F "  cout << \"Calling dtor\" << endl;\n";
    print F "}\n";
    print F "\n";
    print F "\n";



    print F "int $classname\:\:Init(PHCompositeNode *topNode)\n";
    print F "{\n";
    print F "  cout << \"Calling Init\" << endl;\n";
    print F "  return Fun4AllReturnCodes::EVENT_OK;\n";
    print F "}\n";
    print F "\n";

    print F "int $classname\:\:InitRun(PHCompositeNode *topNode)\n";
    print F "{\n";
    print F "  cout << \"Calling InitRun\" << endl;\n";
    print F "  return Fun4AllReturnCodes::EVENT_OK;\n";
    print F "}\n";
    print F "\n";

    print F "int $classname\:\:process_event(PHCompositeNode *topNode)\n";
    print F "{\n";
    print F "  cout << \"Calling process_event\" << endl;\n";
    print F "  return Fun4AllReturnCodes::EVENT_OK;\n";
    print F "}\n";
    print F "\n";

    print F "int $classname\:\:ResetEvent(PHCompositeNode *topNode)\n";
    print F "{\n";
    print F "  cout << \"Resetting internal structures after event processing\" << endl;\n";
    print F "  return Fun4AllReturnCodes::EVENT_OK;\n";
    print F "}\n";
    print F "\n";

    print F "int $classname\:\:EndRun(const int runnumber)\n";
    print F "{\n";
    print F "  cout << \"Ending Run for run \" << runnumber << endl;\n";
    print F "  return Fun4AllReturnCodes::EVENT_OK;\n";
    print F "}\n";
    print F "\n";

    print F "int $classname\:\:End(PHCompositeNode *topNode)\n";
    print F "{\n";
    print F "  cout << \"This is the End...\" << endl;\n";
    print F "  return Fun4AllReturnCodes::EVENT_OK;\n";
    print F "}\n";
    print F "\n";

    print F "int $classname\:\:Reset(PHCompositeNode *topNode)\n";
    print F "{\n";
    print F "  cout << \"being Reset()\" << endl;\n";
    print F "  return Fun4AllReturnCodes::EVENT_OK;\n";
    print F "}\n";
    print F "\n";
    print F "void $classname\:\:Print(const std::string &what) const\n";
    print F "{\n";
    print F "  cout << \"Printing info for \" << what << endl;\n";
    print F "}\n";
    close(F);
}

sub CreateInclude()
{
    my $file = shift;
    open(F,">$file");
    my $includeguard = uc(sprintf("%s_H",$classname));
    print F "// Tell emacs that this is a C++ source\n";
    print F "//  -*- C++ -*-.\n";
    print F "#ifndef $includeguard\n";
    print F "#define $includeguard\n";
    print F "\n";
    print F "#include <fun4all/SubsysReco.h>\n";
    print F "\n";
    print F "#include <string>\n";

    print F "\n";
    print F "class PHCompositeNode;\n";
    print F "\n";
    print F "class $classname : public SubsysReco\n";
    print F "{\n";
    print F " public:\n";
    print F "\n";
    print F "  $classname(const std::string &name = \"$classname\");\n";
    print F "\n";
    print F "  virtual ~$classname();\n";
    print F "\n";
    print F "  /** Called during initialization.\n";
    print F "      Typically this is where you can book histograms, and e.g.\n";
    print F "      register them to Fun4AllServer (so they can be output to file\n";
    print F "      using Fun4AllServer::dumpHistos() method).\n";
    print F "   */\n";
    print F "  int Init(PHCompositeNode *topNode) override;\n";
    print F "\n";
    print F "  /** Called for first event when run number is known.\n";
    print F "      Typically this is where you may want to fetch data from\n";
    print F "      database, because you know the run number.\n";
    print F "   */\n";
    print F "  int InitRun(PHCompositeNode *topNode) override;\n";
    print F "\n";
    print F "  /** Called for each event.\n";
    print F "      This is where you do the real work.\n";
    print F "   */\n";
    print F "  int process_event(PHCompositeNode *topNode) override;\n";
    print F "\n";
    print F "  /// Clean up internals after each event.\n";
    print F "  int ResetEvent(PHCompositeNode *topNode) override;\n";
    print F "\n";
    print F "  /// Called at the end of each run.\n";
    print F "  int EndRun(const int runnumber) override;\n";
    print F "\n";
    print F "  /// Called at the end of all processing.\n";
    print F "  int End(PHCompositeNode *topNode) override;\n";
    print F "\n";

    print F "  /// Reset\n";
    print F "  int Reset(PHCompositeNode * /*topNode*/) override;\n";
    print F "\n";

    print F "  void Print(const std::string &what = \"ALL\") const override;\n";
    print F "\n";

    print F " private:\n";
    print F "};\n";
    print F "\n";

    print F "#endif // $includeguard\n";
    close(F);
}

sub CreateAutogen()
{
    open(F,">autogen.sh");
    print F "#!/bin/sh\n";
    print F "srcdir=`dirname \$0`\n";
    print F "test -z \"\$srcdir\" && srcdir=.\n";
    print F "\n";
    print F "(cd \$srcdir; aclocal -I \${OFFLINE_MAIN}/share;\\\n";
    print F "libtoolize --force; automake -a --add-missing; autoconf)\n";
    print F "\n";
    print F "\$srcdir/configure  \"\$\@\"\n";
    close(F);
    chmod 0755, "autogen.sh";
}

sub CreateConfigure()
{
    my $loname = lc($classname);
    open(F,">configure.ac");
    print F "AC_INIT($loname,[1.00])\n";
    print F "AC_CONFIG_SRCDIR([configure.ac])\n";
    print F "\n";

    print F "AM_INIT_AUTOMAKE\n";
    print F "AC_PROG_CXX(CC g++)\n";
    print F "\n";

    print F "LT_INIT([disable-static])\n";
    print F "\n";

    print F "dnl   no point in suppressing warnings people should \n";
    print F "dnl   at least see them, so here we go for g++: -Wall\n";
    print F "if test \$ac_cv_prog_gxx = yes; then\n";
    print F "   CXXFLAGS=\"\$CXXFLAGS -Wall -Werror\"\n";
    print F "fi\n";
    print F "\n";

    print F "AC_CONFIG_FILES([Makefile])\n";
    print F "AC_OUTPUT\n";
    close(F);
}

sub CreateMakefile()
{
    open(F,">Makefile.am");
    print F "AUTOMAKE_OPTIONS = foreign\n";
    print F "\n";

    print F "AM_CPPFLAGS = \\\n";
    print F "  -I\$(includedir) \\\n";
    print F "  -I\$(OFFLINE_MAIN)/include \\\n";
    print F "  -I\$(ROOTSYS)/include\n";
    print F "\n";

    print F "AM_LDFLAGS = \\\n";
    print F "  -L\$(libdir) \\\n";
    print F "  -L\$(OFFLINE_MAIN)/lib\n";
    print F "\n";

    print F "pkginclude_HEADERS = \\\n";
    print F "  $classname.h\n";
    print F "\n";

    print F "lib_LTLIBRARIES = \\\n";
    print F "  lib$classname.la\n";
    print F "\n";

    print F "lib${classname}_la_SOURCES = \\\n";
    print F "  $classname.cc\n";
    print F "\n";

    print F "lib${classname}_la_LIBADD = \\\n";
    print F "  -lphool \\\n";
    print F "  -lSubsysReco\n";
    print F "\n";

    print F "BUILT_SOURCES = testexternals.cc\n";
    print F "\n";

    print F "noinst_PROGRAMS = \\\n";
    print F "  testexternals\n";
    print F "\n";

    print F "testexternals_SOURCES = testexternals.cc\n";
    print F "testexternals_LDADD   = lib$classname.la\n";
    print F "\n";

    print F "testexternals.cc:\n";
    print F "\techo \"//*** this is a generated file. Do not commit, do not edit\" > \$\@\n";
    print F "\techo \"int main()\" >> \$\@\n";
    print F "\techo \"{\" >> \$\@\n";
    print F "\techo \"  return 0;\" >> \$\@\n";
    print F "\techo \"}\" >> \$\@\n";
    print F "\n";

    print F "clean-local:\n";
    print F "\trm -f \$(BUILT_SOURCES)\n";
    close(F);
}

