#!/usr/bin/env python3
import argparse
import json
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

BASEDIR = Path(__file__).parent.parent
DEFAULT_OUT = str(BASEDIR / "data" / "access_token.json")
DEFAULT_SSO = str(BASEDIR / "data" / "input_sso.json")


def fetch_access_token(auth_base_url, client_id, client_secret, scope, timeout):
    token_url = f"{auth_base_url.rstrip('/')}/protocol/openid-connect/token"
    payload = {
        "grant_type": "client_credentials",
        "client_id": client_id,
        "client_secret": client_secret,
        "scope": scope,
    }

    body = urlencode(payload).encode("utf-8")
    request = Request(
        token_url,
        data=body,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )

    try:
        with urlopen(request, timeout=timeout) as response:
            response_text = response.read().decode("utf-8", errors="replace")
            status_code = response.getcode()
    except HTTPError as exc:
        error_body = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"Token request failed ({exc.code}): {error_body}") from exc
    except URLError as exc:
        raise RuntimeError(f"Token request failed: {exc.reason}") from exc

    if status_code < 200 or status_code >= 300:
        raise RuntimeError(f"Token request failed ({status_code}): {response_text}")

    try:
        token_data = json.loads(response_text)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"Token response is not valid JSON: {response_text}") from exc

    access_token = token_data.get("access_token")
    if not access_token:
        raise RuntimeError(f"No access_token found in response: {response_text}")

    issued_at = datetime.now(timezone.utc)
    expires_in = token_data.get("expires_in")
    expires_at = None
    if isinstance(expires_in, int):
        expires_at = (issued_at + timedelta(seconds=expires_in)).isoformat()

    return {
        "token_url": token_url,
        "issued_at": issued_at.isoformat(),
        "expires_at": expires_at,
        "access_token": access_token,
        "token_type": token_data.get("token_type", "Bearer"),
        "expires_in": expires_in,
        "scope": token_data.get("scope", scope),
        "raw_response": token_data,
    }


def main():
    parser = argparse.ArgumentParser(description="Get OAuth2 access token using client credentials")
    parser.add_argument("--sso", default=DEFAULT_SSO, help=f"Path to SSO input JSON (default: {DEFAULT_SSO})")
    parser.add_argument("--scope", default="openid", help="OAuth scope (default: openid)")
    parser.add_argument("--timeout", type=int, default=30, help="Request timeout in seconds (default: 30)")
    parser.add_argument(
        "--format",
        choices=["token", "bearer", "json"],
        default=None,
        help="Optional CLI output format: token, bearer, or json (default: no CLI output)",
    )
    parser.add_argument("--out", default=DEFAULT_OUT, help=f"Path to write full JSON output (default: {DEFAULT_OUT})")
    args = parser.parse_args()

    sso_path = Path(args.sso)
    if not sso_path.exists():
        print(f"ERROR: SSO input file not found: {sso_path}", file=sys.stderr)
        sys.exit(1)
    with open(sso_path, "r") as sso_handle:
        sso = json.load(sso_handle)
    for key in ("auth_base_url", "client_id", "client_secret"):
        if not sso.get(key):
            print(f"ERROR: Missing '{key}' in {sso_path}", file=sys.stderr)
            sys.exit(1)

    try:
        token_result = fetch_access_token(
            auth_base_url=sso["auth_base_url"],
            client_id=sso["client_id"],
            client_secret=sso["client_secret"],
            scope=args.scope,
            timeout=args.timeout,
        )
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as output_handle:
            json.dump(token_result, output_handle, indent=2)

    if args.format:
        if args.format == "token":
            print(token_result["access_token"])
        elif args.format == "bearer":
            print(f"Bearer {token_result['access_token']}")
        elif args.format == "json":
            print(json.dumps(token_result, indent=2))


if __name__ == "__main__":
    main()
