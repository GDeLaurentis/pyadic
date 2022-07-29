#!/usr/bin/env python
import re

repo_owner = "GDeLaurentis"
repo_name = "pyadic"


def update_coverage_badge():
    """Updates the coverage badge in README.md"""
    # read value from coverage.txt
    with open('coverage.txt', 'r') as f:
        content = f.read()
    content = content.split("TOTAL   ")[1].split("\n\n")[0]
    coverage_percentage = int(re.findall(r"(\d+)%", content)[0])
    # obtain colour of badge
    coverage_steps = {
        100: "brightgreen",
        90: "green",
        80: "greenyellow",
        70: "yellow",
        50: "orange",
        0: "red",
    }
    for key in sorted(coverage_steps.keys())[::-1]:
        if coverage_percentage >= key:
            break
    coverage_colour = coverage_steps[key]
    # update README.md
    with open('README.md', 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if "![Coverage]" in line:
            lines[i] = '[![Coverage](https://img.shields.io/badge/Coverage-{}%25-{}?labelColor=2a2f35)](https://github.com/{}/{}/actions)\n'.format(
                coverage_percentage, coverage_colour, repo_owner, repo_name)
    with open('README.md', 'w') as f:
        f.write("".join(lines))


if __name__ == "__main__":
    update_coverage_badge()
